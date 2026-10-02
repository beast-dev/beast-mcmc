package test.dr.inference.model;

import dr.evomodel.treedatalikelihood.continuous.canonical.gradient.CanonicalSelectionGradientProjector;
import dr.inference.model.BlockDiagonalPolarStableMatrixParameter;
import dr.inference.model.CompoundParameter;
import dr.inference.model.GivensRotationMatrixParameter;
import dr.inference.model.InvertibleMatrixParametrization;
import dr.inference.model.OrthogonalMatrixProvider;
import dr.inference.model.Parameter;
import dr.inference.model.QRUnitUpperMatrixParametrization;
import test.dr.math.MathTestCase;

public class QRUnitUpperMatrixParametrizationTest extends MathTestCase {

    public void testZeroParametersGiveIdentityAndExposeNativeState() {
        final Parameter angles = new Parameter.Default("angles", new double[3]);
        final Parameter upper = new Parameter.Default("upper", new double[3]);
        final GivensRotationMatrixParameter orthogonal =
                new GivensRotationMatrixParameter("orthogonal", 3, angles);
        final QRUnitUpperMatrixParametrization matrix =
                new QRUnitUpperMatrixParametrization("qr", orthogonal, upper);

        assertTrue(matrix instanceof InvertibleMatrixParametrization);
        assertTrue(matrix.isInvertible());
        assertTrue(matrix.supportsNativeParameter(angles));
        assertTrue(matrix.supportsNativeParameter(upper));
        assertTrue(matrix.supportsNativeParameter(matrix.getNativeParameter()));
        assertTrue(matrix.getOrthogonalProvider() == orthogonal);
        assertTrue(matrix.getOrthogonalMatrixParameter() == orthogonal);
        assertTrue(matrix.getOrthogonalParameter() == angles);
        assertTrue(matrix.getUpperParameter() == upper);
        assertTrue(matrix.getCheckpointableStateParameter() == matrix.getNativeParameter());
        assertEquals(6, matrix.getNativeDimension());
        assertEquals(0.0, matrix.getLogAbsDeterminant(), 1e-15);
        assertEquals(1.0, matrix.getDeterminantSign(), 1e-15);
        assertEquals(1.0, matrix.getDeterminant(), 1e-15);

        for (int row = 0; row < 3; row++) {
            for (int col = 0; col < 3; col++) {
                assertEquals(row == col ? 1.0 : 0.0, matrix.getParameterValue(row, col), 1e-15);
            }
        }
    }

    public void testMatrixEqualsOrthogonalTimesUnitUpper() {
        final Parameter angles = new Parameter.Default("angles", new double[]{0.2, -0.4, 0.7});
        final Parameter upper = new Parameter.Default("upper", new double[]{-0.1, 0.6, 0.25});
        final GivensRotationMatrixParameter orthogonal =
                new GivensRotationMatrixParameter("orthogonal", 3, angles);
        final QRUnitUpperMatrixParametrization matrix =
                new QRUnitUpperMatrixParametrization("qr", orthogonal, upper);
        final double[] q = new double[9];
        final double[] unitUpper = new double[]{
                1.0, -0.1, 0.6,
                0.0, 1.0, 0.25,
                0.0, 0.0, 1.0
        };

        orthogonal.fillOrthogonalMatrix(q);

        for (int row = 0; row < 3; row++) {
            for (int col = 0; col < 3; col++) {
                double expected = 0.0;
                for (int k = 0; k < 3; k++) {
                    expected += q[row * 3 + k] * unitUpper[k * 3 + col];
                }
                assertEquals(expected, matrix.getParameterValue(row, col), 1e-15);
            }
        }
    }

    public void testInverseMultipliesToIdentity() {
        final Parameter angles = new Parameter.Default("angles", new double[]{0.2, -0.4, 0.7});
        final Parameter upper = new Parameter.Default("upper", new double[]{-0.1, 0.6, 0.25});
        final QRUnitUpperMatrixParametrization matrix =
                new QRUnitUpperMatrixParametrization(
                        "qr",
                        new GivensRotationMatrixParameter("orthogonal", 3, angles),
                        upper);

        final double[] inverse = new double[9];
        matrix.fillInverse(inverse);

        assertProductIsIdentity(matrix, inverse);
    }

    public void testPullBackGradientMatchesFiniteDifference() {
        final Parameter angles = new Parameter.Default("angles", new double[]{0.2, -0.4, 0.7});
        final Parameter upper = new Parameter.Default("upper", new double[]{-0.1, 0.6, 0.25});
        final QRUnitUpperMatrixParametrization matrix =
                new QRUnitUpperMatrixParametrization(
                        "qr",
                        new GivensRotationMatrixParameter("orthogonal", 3, angles),
                        upper);
        final double[] gradientR = new double[]{
                0.7, -1.1, 0.4,
                0.3, 0.9, -0.8,
                -0.2, 1.4, 0.6
        };

        final double[] analytic = new double[6];
        matrix.fillPullBackGradientFlat(gradientR, 3, analytic);

        final Parameter[] parameters = new Parameter[]{angles, upper};
        int offset = 0;
        for (Parameter parameter : parameters) {
            for (int i = 0; i < parameter.getDimension(); i++) {
                assertEquals(
                        "finite difference gradient for " + parameter.getId() + "[" + i + "]",
                        finiteDifference(matrix, parameter, i, gradientR),
                        analytic[offset + i],
                        1e-6);
            }
            offset += parameter.getDimension();
        }
    }

    public void testPullBackGradientForNativeSubParameters() {
        final Parameter angles = new Parameter.Default("angles", new double[]{0.2, -0.4, 0.7});
        final Parameter upper = new Parameter.Default("upper", new double[]{-0.1, 0.6, 0.25});
        final QRUnitUpperMatrixParametrization matrix =
                new QRUnitUpperMatrixParametrization(
                        "qr",
                        new GivensRotationMatrixParameter("orthogonal", 3, angles),
                        upper);
        final double[] gradientR = new double[]{
                0.7, -1.1, 0.4,
                0.3, 0.9, -0.8,
                -0.2, 1.4, 0.6
        };

        final double[] full = new double[6];
        matrix.fillPullBackGradientFlat(gradientR, 3, full);

        final double[] angleGradient = matrix.pullBackGradientForParameter(angles, gradientR, 3);
        final double[] upperGradient = matrix.pullBackGradientForParameter(upper, gradientR, 3);
        assertEquals(3, angleGradient.length);
        assertEquals(3, upperGradient.length);
        for (int i = 0; i < 3; i++) {
            assertEquals(full[i], angleGradient[i], 1e-15);
            assertEquals(full[3 + i], upperGradient[i], 1e-15);
        }

        final CompoundParameter nativeCompound = new CompoundParameter("qrNative");
        nativeCompound.addParameter(angles);
        nativeCompound.addParameter(upper);
        final double[] compoundGradient = matrix.pullBackGradientForParameter(nativeCompound, gradientR, 3);
        for (int i = 0; i < full.length; i++) {
            assertEquals(full[i], compoundGradient[i], 1e-15);
        }
    }

    public void testCanonicalProjectorRoutesGradientToQrNativeParameters() {
        final Parameter angles = new Parameter.Default("angles", new double[]{0.2});
        final Parameter upper = new Parameter.Default("upper", new double[]{-0.1});
        final QRUnitUpperMatrixParametrization rotation =
                new QRUnitUpperMatrixParametrization(
                        "qr",
                        new GivensRotationMatrixParameter("orthogonal", 2, angles),
                        upper);
        final BlockDiagonalPolarStableMatrixParameter block =
                new BlockDiagonalPolarStableMatrixParameter(
                        "block",
                        rotation,
                        new Parameter.Default(0),
                        new Parameter.Default("rho", new double[]{0.9}),
                        new Parameter.Default("theta", new double[]{0.3}),
                        new Parameter.Default("t", new double[]{0.1}));
        final CompoundParameter qrNative = new CompoundParameter("qrNative");
        qrNative.addParameter(angles);
        qrNative.addParameter(upper);
        final double[] nativeBlockGradient = new double[]{1.0, 2.0, 3.0};
        final double[] gradientR = new double[]{
                0.7, -1.1,
                0.3, 0.9
        };
        final double[] fullRotationGradient = new double[2];
        rotation.fillPullBackGradientFlat(gradientR, 2, fullRotationGradient);

        final CompoundParameter checkpointableState = (CompoundParameter) block.getCheckpointableStateParameter();
        assertTrue(checkpointableState.getParameter(4) == rotation.getNativeParameter());

        final double[] angleGradient = CanonicalSelectionGradientProjector.assembleBlockGradientResultFlat(
                2, angles, block, nativeBlockGradient, gradientR);
        final double[] upperGradient = CanonicalSelectionGradientProjector.assembleBlockGradientResultFlat(
                2, upper, block, nativeBlockGradient, gradientR);
        final double[] compoundGradient = CanonicalSelectionGradientProjector.assembleBlockGradientResultFlat(
                2, qrNative, block, nativeBlockGradient, gradientR);

        assertEquals(1, angleGradient.length);
        assertEquals(1, upperGradient.length);
        assertEquals(2, compoundGradient.length);
        assertEquals(fullRotationGradient[0], angleGradient[0], 1e-15);
        assertEquals(fullRotationGradient[1], upperGradient[0], 1e-15);
        for (int i = 0; i < fullRotationGradient.length; i++) {
            assertEquals(fullRotationGradient[i], compoundGradient[i], 1e-15);
        }
    }

    public void testDimensionChecksRejectWrongUpperSizeAndRawMatrix() {
        try {
            new QRUnitUpperMatrixParametrization(
                    "badUpper",
                    new GivensRotationMatrixParameter(
                            "orthogonal",
                            3,
                            new Parameter.Default("angles", new double[3])),
                    new Parameter.Default("upper", new double[2]));
            fail("Expected IllegalArgumentException for upper dimension mismatch");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage().contains("unitUpper"));
        }

        try {
            new QRUnitUpperMatrixParametrization(
                    "badOrthogonal",
                    new dr.inference.model.MatrixParameter("raw", 3, 3),
                    new Parameter.Default("upper", new double[3]));
            fail("Expected IllegalArgumentException for non-orthogonal matrix");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage().contains("OrthogonalMatrixProvider"));
        }
    }

    private static void assertProductIsIdentity(final QRUnitUpperMatrixParametrization matrix,
                                                final double[] inverse) {
        final int dim = matrix.getRowDimension();
        for (int row = 0; row < dim; row++) {
            for (int col = 0; col < dim; col++) {
                double product = 0.0;
                for (int k = 0; k < dim; k++) {
                    product += matrix.getParameterValue(row, k) * inverse[k * dim + col];
                }
                assertEquals(row == col ? 1.0 : 0.0, product, 1e-14);
            }
        }
    }

    private static double finiteDifference(final QRUnitUpperMatrixParametrization matrix,
                                           final Parameter parameter,
                                           final int index,
                                           final double[] gradientR) {
        final double h = 1e-7;
        final double original = parameter.getParameterValue(index);
        parameter.setParameterValue(index, original + h);
        final double plus = objective(matrix, gradientR);
        parameter.setParameterValue(index, original - h);
        final double minus = objective(matrix, gradientR);
        parameter.setParameterValue(index, original);
        return (plus - minus) / (2.0 * h);
    }

    private static double objective(final QRUnitUpperMatrixParametrization matrix,
                                    final double[] gradientR) {
        double objective = 0.0;
        final int dim = matrix.getRowDimension();
        for (int row = 0; row < dim; row++) {
            for (int col = 0; col < dim; col++) {
                objective += gradientR[row * dim + col] * matrix.getParameterValue(row, col);
            }
        }
        return objective;
    }
}
