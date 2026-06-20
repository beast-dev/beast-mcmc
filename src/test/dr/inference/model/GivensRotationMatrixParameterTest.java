package test.dr.inference.model;

import dr.evomodel.treedatalikelihood.continuous.canonical.gradient.CanonicalSelectionGradientProjector;
import dr.inference.model.BlockDiagonalPolarStableMatrixParameter;
import dr.inference.model.GivensRotationMatrixParameter;
import dr.inference.model.InvertibleMatrixParametrization;
import dr.inference.model.OrthogonalMatrixProvider;
import dr.inference.model.Parameter;
import test.dr.math.MathTestCase;

public class GivensRotationMatrixParameterTest extends MathTestCase {

    public void testGivensRotationIsInvertibleParametrization() {
        final Parameter angles = new Parameter.Default("angles", new double[]{0.2, -0.4, 0.7});
        final GivensRotationMatrixParameter rotation =
                new GivensRotationMatrixParameter("rotation", 3, angles);

        assertTrue(rotation instanceof OrthogonalMatrixProvider);
        assertTrue(rotation instanceof InvertibleMatrixParametrization);
        assertTrue(rotation.isInvertible());
        assertTrue(rotation.supportsNativeParameter(angles));
        assertTrue(rotation.supportsNativeParameter(rotation.getNativeParameter()));
        assertTrue(rotation.getNativeParameter().getParameter(0) == angles);
        assertEquals(3, rotation.getNativeDimension());
        assertEquals(0.0, rotation.getLogAbsDeterminant(), 1e-15);
        assertEquals(1.0, rotation.getDeterminantSign(), 1e-15);
        assertEquals(1.0, rotation.getDeterminant(), 1e-15);
    }

    public void testInverseIsOrthogonalTranspose() {
        final Parameter angles = new Parameter.Default("angles", new double[]{0.2, -0.4, 0.7});
        final GivensRotationMatrixParameter rotation =
                new GivensRotationMatrixParameter("rotation", 3, angles);

        final double[] inverse = new double[9];
        final double[] transpose = new double[9];
        rotation.fillInverse(inverse);
        rotation.fillOrthogonalTranspose(transpose);

        for (int i = 0; i < inverse.length; i++) {
            assertEquals(transpose[i], inverse[i], 1e-15);
        }
        assertProductIsIdentity(rotation, inverse);
    }

    public void testNativeGradientMatchesOrthogonalAngleGradient() {
        final Parameter angles = new Parameter.Default("angles", new double[]{0.2, -0.4, 0.7});
        final GivensRotationMatrixParameter rotation =
                new GivensRotationMatrixParameter("rotation", 3, angles);
        final double[] gradientR = new double[]{
                0.7, -1.1, 0.4,
                0.3, 0.9, -0.8,
                -0.2, 1.4, 0.6
        };

        final double[] orthogonalGradient = new double[3];
        rotation.fillPullBackGradientFlat(gradientR, 3, orthogonalGradient);

        final double[] nativeGradient = new double[3];
        rotation.fillPullBackGradientForParameter(angles, gradientR, 3, nativeGradient, 0);

        final double[] nativeCompoundGradient = new double[3];
        rotation.fillPullBackGradientForParameter(
                rotation.getNativeParameter(), gradientR, 3, nativeCompoundGradient, 0);

        for (int i = 0; i < orthogonalGradient.length; i++) {
            assertEquals(orthogonalGradient[i], nativeGradient[i], 1e-15);
            assertEquals(orthogonalGradient[i], nativeCompoundGradient[i], 1e-15);
            assertEquals(
                    "finite difference gradient for angles[" + i + "]",
                    finiteDifference(rotation, angles, i, gradientR),
                    nativeGradient[i],
                    1e-6);
        }
    }

    public void testCanonicalProjectorRoutesGivensNativeCompound() {
        final Parameter angles = new Parameter.Default("angles", new double[]{0.2});
        final GivensRotationMatrixParameter rotation =
                new GivensRotationMatrixParameter("rotation", 2, angles);
        final BlockDiagonalPolarStableMatrixParameter block =
                new BlockDiagonalPolarStableMatrixParameter(
                        "block",
                        rotation,
                        new Parameter.Default(0),
                        new Parameter.Default("rho", new double[]{0.9}),
                        new Parameter.Default("theta", new double[]{0.3}),
                        new Parameter.Default("t", new double[]{0.1}));
        final double[] nativeBlockGradient = new double[]{1.0, 2.0, 3.0};
        final double[] gradientR = new double[]{
                0.7, -1.1,
                0.3, 0.9
        };
        final double[] expected = new double[1];
        rotation.fillPullBackGradientForParameter(rotation.getNativeParameter(), gradientR, 2, expected, 0);

        final double[] actual = CanonicalSelectionGradientProjector.assembleBlockGradientResultFlat(
                2, rotation.getNativeParameter(), block, nativeBlockGradient, gradientR);

        assertEquals(1, actual.length);
        assertEquals(expected[0], actual[0], 1e-15);
    }

    private static void assertProductIsIdentity(final GivensRotationMatrixParameter rotation,
                                                final double[] inverse) {
        final int dim = rotation.getRowDimension();
        for (int row = 0; row < dim; row++) {
            for (int col = 0; col < dim; col++) {
                double product = 0.0;
                for (int k = 0; k < dim; k++) {
                    product += rotation.getParameterValue(row, k) * inverse[k * dim + col];
                }
                assertEquals(row == col ? 1.0 : 0.0, product, 1e-14);
            }
        }
    }

    private static double finiteDifference(final GivensRotationMatrixParameter rotation,
                                           final Parameter angles,
                                           final int index,
                                           final double[] gradientR) {
        final double h = 1e-7;
        final double original = angles.getParameterValue(index);
        angles.setParameterValue(index, original + h);
        final double plus = objective(rotation, gradientR);
        angles.setParameterValue(index, original - h);
        final double minus = objective(rotation, gradientR);
        angles.setParameterValue(index, original);
        return (plus - minus) / (2.0 * h);
    }

    private static double objective(final GivensRotationMatrixParameter rotation,
                                    final double[] gradientR) {
        double objective = 0.0;
        final int dim = rotation.getRowDimension();
        for (int row = 0; row < dim; row++) {
            for (int col = 0; col < dim; col++) {
                objective += gradientR[row * dim + col] * rotation.getParameterValue(row, col);
            }
        }
        return objective;
    }
}
