package test.dr.inference.model;

import dr.evomodel.treedatalikelihood.continuous.canonical.gradient.CanonicalSelectionGradientProjector;
import dr.inference.model.BlockDiagonalPolarStableMatrixParameter;
import dr.inference.model.CompoundParameter;
import dr.inference.model.InvertibleMatrixParametrization;
import dr.inference.model.LDUInvertibleMatrixParametrization;
import dr.inference.model.Parameter;
import test.dr.math.MathTestCase;

import java.util.Collections;

public class LDUInvertibleMatrixParametrizationTest extends MathTestCase {

    public void testZeroParametersGiveIdentity() {
        final LDUInvertibleMatrixParametrization matrix =
                new LDUInvertibleMatrixParametrization(
                        "ldu",
                        3,
                        new Parameter.Default("lower", new double[3]),
                        new Parameter.Default("logDiagonal", new double[3]),
                        new Parameter.Default("upper", new double[3]));

        for (int row = 0; row < 3; row++) {
            for (int col = 0; col < 3; col++) {
                assertEquals(row == col ? 1.0 : 0.0, matrix.getParameterValue(row, col), 1e-15);
            }
        }
        assertEquals(0.0, matrix.getLogAbsDeterminant(), 1e-15);
        assertEquals(1.0, matrix.getDeterminantSign(), 1e-15);
        assertEquals(1.0, matrix.getDeterminant(), 1e-15);
    }

    public void testExpectedThreeByThreeMatrixWithFixedDiagonalSigns() {
        final LDUInvertibleMatrixParametrization matrix =
                new LDUInvertibleMatrixParametrization(
                        "ldu",
                        3,
                        new Parameter.Default("lower", new double[]{0.2, -0.3, 0.4}),
                        new Parameter.Default("logDiagonal", new double[]{
                                Math.log(2.0), Math.log(3.0), Math.log(5.0)}),
                        new Parameter.Default("upper", new double[]{-0.1, 0.6, 0.25}),
                        new double[]{1.0, -1.0, 1.0});

        assertEquals(2.0, matrix.getParameterValue(0, 0), 1e-15);
        assertEquals(-0.2, matrix.getParameterValue(0, 1), 1e-15);
        assertEquals(1.2, matrix.getParameterValue(0, 2), 1e-15);
        assertEquals(0.4, matrix.getParameterValue(1, 0), 1e-15);
        assertEquals(-3.04, matrix.getParameterValue(1, 1), 1e-15);
        assertEquals(-0.51, matrix.getParameterValue(1, 2), 1e-15);
        assertEquals(-0.6, matrix.getParameterValue(2, 0), 1e-15);
        assertEquals(-1.14, matrix.getParameterValue(2, 1), 1e-15);
        assertEquals(4.34, matrix.getParameterValue(2, 2), 1e-15);

        assertEquals(Math.log(30.0), matrix.getLogAbsDeterminant(), 1e-15);
        assertEquals(-1.0, matrix.getDeterminantSign(), 1e-15);
        assertEquals(-30.0, matrix.getDeterminant(), 1e-12);
    }

    public void testInverseMultipliesToIdentity() {
        final LDUInvertibleMatrixParametrization matrix =
                new LDUInvertibleMatrixParametrization(
                        "ldu",
                        3,
                        new Parameter.Default("lower", new double[]{0.2, -0.3, 0.4}),
                        new Parameter.Default("logDiagonal", new double[]{
                                Math.log(2.0), Math.log(3.0), Math.log(5.0)}),
                        new Parameter.Default("upper", new double[]{-0.1, 0.6, 0.25}));

        final double[] inverse = new double[9];
        matrix.fillInverse(inverse);

        for (int row = 0; row < 3; row++) {
            for (int col = 0; col < 3; col++) {
                double product = 0.0;
                for (int k = 0; k < 3; k++) {
                    product += matrix.getParameterValue(row, k) * inverse[k * 3 + col];
                }
                assertEquals(row == col ? 1.0 : 0.0, product, 1e-14);
            }
        }
    }

    public void testPullBackGradientMatchesFiniteDifference() {
        final Parameter lower = new Parameter.Default("lower", new double[]{0.2, -0.3, 0.4});
        final Parameter logDiagonal = new Parameter.Default("logDiagonal", new double[]{
                Math.log(2.0), Math.log(3.0), Math.log(5.0)});
        final Parameter upper = new Parameter.Default("upper", new double[]{-0.1, 0.6, 0.25});
        final LDUInvertibleMatrixParametrization matrix =
                new LDUInvertibleMatrixParametrization(
                        "ldu", 3, lower, logDiagonal, upper);
        final double[] gradientR = new double[]{
                0.7, -1.1, 0.4,
                0.3, 0.9, -0.8,
                -0.2, 1.4, 0.6
        };

        final double[] analytic = new double[9];
        matrix.fillPullBackGradientFlat(gradientR, 3, analytic);

        final Parameter[] parameters = new Parameter[]{lower, logDiagonal, upper};
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
        final Parameter lower = new Parameter.Default("lower", new double[]{0.2, -0.3, 0.4});
        final Parameter logDiagonal = new Parameter.Default("logDiagonal", new double[]{
                Math.log(2.0), Math.log(3.0), Math.log(5.0)});
        final Parameter upper = new Parameter.Default("upper", new double[]{-0.1, 0.6, 0.25});
        final LDUInvertibleMatrixParametrization matrix =
                new LDUInvertibleMatrixParametrization(
                        "ldu", 3, lower, logDiagonal, upper);
        final double[] gradientR = new double[]{
                0.7, -1.1, 0.4,
                0.3, 0.9, -0.8,
                -0.2, 1.4, 0.6
        };

        final double[] full = new double[9];
        matrix.fillPullBackGradientFlat(gradientR, 3, full);

        assertTrue(matrix.supportsNativeParameter(matrix.getNativeParameter()));
        assertTrue(matrix.getNativeParameter().getParameter(0) == lower);
        assertTrue(matrix.getNativeParameter().getParameter(1) == logDiagonal);
        assertTrue(matrix.getNativeParameter().getParameter(2) == upper);

        final double[] lowerGradient = matrix.pullBackGradientForParameter(lower, gradientR, 3);
        final double[] logDiagonalGradient = matrix.pullBackGradientForParameter(logDiagonal, gradientR, 3);
        final double[] upperGradient = matrix.pullBackGradientForParameter(upper, gradientR, 3);
        final double[] lowerGradientFilled = new double[matrix.getPullBackGradientDimension(lower)];
        matrix.fillPullBackGradientForParameter(lower, gradientR, 3, lowerGradientFilled, 0);

        assertEquals(3, lowerGradient.length);
        assertEquals(3, logDiagonalGradient.length);
        assertEquals(3, upperGradient.length);
        assertEquals(3, lowerGradientFilled.length);
        for (int i = 0; i < 3; i++) {
            assertEquals(full[i], lowerGradient[i], 1e-15);
            assertEquals(full[i], lowerGradientFilled[i], 1e-15);
            assertEquals(full[3 + i], logDiagonalGradient[i], 1e-15);
            assertEquals(full[6 + i], upperGradient[i], 1e-15);
        }

        final CompoundParameter nativeCompound = new CompoundParameter("native");
        nativeCompound.addParameter(lower);
        nativeCompound.addParameter(logDiagonal);
        nativeCompound.addParameter(upper);
        final double[] compoundGradient = matrix.pullBackGradientForParameter(nativeCompound, gradientR, 3);
        final double[] compoundGradientFilled = new double[matrix.getPullBackGradientDimension(nativeCompound)];
        matrix.fillPullBackGradientForParameter(nativeCompound, gradientR, 3, compoundGradientFilled, 0);
        for (int i = 0; i < full.length; i++) {
            assertEquals(full[i], compoundGradient[i], 1e-15);
            assertEquals(full[i], compoundGradientFilled[i], 1e-15);
        }

        final double[] internalCompoundGradient =
                matrix.pullBackGradientForParameter(matrix.getNativeParameter(), gradientR, 3);
        for (int i = 0; i < full.length; i++) {
            assertEquals(full[i], internalCompoundGradient[i], 1e-15);
        }
    }

    public void testCanonicalProjectorRoutesGradientToLduNativeParameters() {
        final Parameter lower = new Parameter.Default("lower", new double[]{0.2});
        final Parameter logDiagonal = new Parameter.Default("logDiagonal", new double[]{
                Math.log(2.0), Math.log(3.0)});
        final Parameter upper = new Parameter.Default("upper", new double[]{-0.1});
        final LDUInvertibleMatrixParametrization rotation =
                new LDUInvertibleMatrixParametrization(
                        "ldu", 2, lower, logDiagonal, upper);
        final BlockDiagonalPolarStableMatrixParameter block =
                new BlockDiagonalPolarStableMatrixParameter(
                        "block",
                        rotation,
                        new Parameter.Default(0),
                        new Parameter.Default("rho", new double[]{0.9}),
                        new Parameter.Default("theta", new double[]{0.3}),
                        new Parameter.Default("t", new double[]{0.1}));
        final CompoundParameter lduNative = new CompoundParameter("lduNative");
        lduNative.addParameter(lower);
        lduNative.addParameter(logDiagonal);
        lduNative.addParameter(upper);
        final double[] nativeBlockGradient = new double[]{1.0, 2.0, 3.0};
        final double[] gradientR = new double[]{
                0.7, -1.1,
                0.3, 0.9
        };
        final double[] fullRotationGradient = new double[4];
        rotation.fillPullBackGradientFlat(gradientR, 2, fullRotationGradient);

        final double[] lowerGradient = CanonicalSelectionGradientProjector.assembleBlockGradientResultFlat(
                2, lower, block, nativeBlockGradient, gradientR);
        final double[] logDiagonalGradient = CanonicalSelectionGradientProjector.assembleBlockGradientResultFlat(
                2, logDiagonal, block, nativeBlockGradient, gradientR);
        final double[] upperGradient = CanonicalSelectionGradientProjector.assembleBlockGradientResultFlat(
                2, upper, block, nativeBlockGradient, gradientR);
        final double[] compoundGradient = CanonicalSelectionGradientProjector.assembleBlockGradientResultFlat(
                2, lduNative, block, nativeBlockGradient, gradientR);

        assertEquals(1, lowerGradient.length);
        assertEquals(2, logDiagonalGradient.length);
        assertEquals(1, upperGradient.length);
        assertEquals(4, compoundGradient.length);
        assertEquals(fullRotationGradient[0], lowerGradient[0], 1e-15);
        assertEquals(fullRotationGradient[1], logDiagonalGradient[0], 1e-15);
        assertEquals(fullRotationGradient[2], logDiagonalGradient[1], 1e-15);
        assertEquals(fullRotationGradient[3], upperGradient[0], 1e-15);
        for (int i = 0; i < fullRotationGradient.length; i++) {
            assertEquals(fullRotationGradient[i], compoundGradient[i], 1e-15);
        }
    }

    public void testBlockDiagonalUsesInvertibleMatrixInversePath() {
        final CountingInvertibleMatrixParametrization rotation =
                new CountingInvertibleMatrixParametrization(
                        "countingRotation",
                        new Parameter.Default("native", new double[]{0.0}));
        final BlockDiagonalPolarStableMatrixParameter block =
                new BlockDiagonalPolarStableMatrixParameter(
                        "block",
                        rotation,
                        new Parameter.Default(0),
                        new Parameter.Default("rho", new double[]{0.9}),
                        new Parameter.Default("theta", new double[]{0.3}),
                        new Parameter.Default("t", new double[]{0.1}));

        final double[] r = new double[4];
        final double[] rinv = new double[4];
        block.fillRAndRinv(r, rinv);

        assertEquals(1, rotation.getFillInverseCount());
        assertEquals(2.0, r[0], 1e-15);
        assertEquals(0.0, r[1], 1e-15);
        assertEquals(0.0, r[2], 1e-15);
        assertEquals(0.5, r[3], 1e-15);
        assertEquals(0.5, rinv[0], 1e-15);
        assertEquals(0.0, rinv[1], 1e-15);
        assertEquals(0.0, rinv[2], 1e-15);
        assertEquals(2.0, rinv[3], 1e-15);
    }

    public void testDimensionChecksRejectWrongParameterSizes() {
        try {
            new LDUInvertibleMatrixParametrization(
                    "bad",
                    3,
                    new Parameter.Default("lower", new double[2]),
                    new Parameter.Default("logDiagonal", new double[3]),
                    new Parameter.Default("upper", new double[3]));
            fail("Expected IllegalArgumentException for lower dimension mismatch");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage().contains("lower"));
        }
    }

    private static double finiteDifference(final LDUInvertibleMatrixParametrization matrix,
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

    private static double objective(final LDUInvertibleMatrixParametrization matrix,
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

    private static final class CountingInvertibleMatrixParametrization
            extends InvertibleMatrixParametrization {

        private static final int DIMENSION = 2;
        private final double[] matrix = new double[]{2.0, 0.0, 0.0, 0.5};
        private int fillInverseCount = 0;

        private CountingInvertibleMatrixParametrization(final String name,
                                                       final Parameter nativeParameter) {
            super(name, DIMENSION, Collections.singletonList(nativeParameter));
        }

        @Override
        protected double computeEntry(final int row, final int col) {
            return matrix[row * DIMENSION + col];
        }

        @Override
        protected double getCachedValue(final int row, final int col) {
            return matrix[row * DIMENSION + col];
        }

        @Override
        protected void updateCache() {
            // The test matrix is fixed; there is no derived cache to refresh.
        }

        @Override
        public double getLogAbsDeterminant() {
            return 0.0;
        }

        @Override
        public double getDeterminantSign() {
            return 1.0;
        }

        @Override
        public void fillInverse(final double[] inverseRowMajor) {
            fillInverseCount++;
            inverseRowMajor[0] = 0.5;
            inverseRowMajor[1] = 0.0;
            inverseRowMajor[2] = 0.0;
            inverseRowMajor[3] = 2.0;
        }

        @Override
        public void fillPullBackGradientFlat(final double[] gradientWrtMatrixRowMajor,
                                             final int dimension,
                                             final double[] out,
                                             final int offset) {
            out[offset] = 0.0;
        }

        private int getFillInverseCount() {
            return fillInverseCount;
        }
    }
}
