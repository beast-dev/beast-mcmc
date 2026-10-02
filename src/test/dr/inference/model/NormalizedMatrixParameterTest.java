package test.dr.inference.model;

import dr.evomodel.treedatalikelihood.continuous.canonical.gradient.CanonicalSelectionGradientProjector;
import dr.inference.model.BlockDiagonalPolarStableMatrixParameter;
import dr.inference.model.MatrixParameter;
import dr.inference.model.NativeMatrixGradientProvider;
import dr.inference.model.NormalizedMatrixParameter;
import dr.inference.model.Parameter;
import test.dr.math.MathTestCase;

public class NormalizedMatrixParameterTest extends MathTestCase {

    public void testGlobalNormalizationHasUnitFrobeniusNorm() {
        final MatrixParameter raw = makeMatrixParameter("raw", new double[][]{
                {1.0, 2.0},
                {3.0, 4.0}
        });
        final NormalizedMatrixParameter matrix = new NormalizedMatrixParameter(
                "normalized", raw, NormalizedMatrixParameter.Normalization.GLOBAL);

        assertTrue(matrix instanceof NativeMatrixGradientProvider);
        assertTrue(matrix.supportsNativeParameter(raw));
        assertTrue(matrix.getNativeParameter() == raw);
        assertEquals(Math.sqrt(30.0), matrix.getGlobalNorm(), 1e-15);
        assertEquals(1.0, frobeniusNorm(matrix), 1e-15);
    }

    public void testColumnNormalizationHasUnitColumns() {
        final MatrixParameter raw = makeMatrixParameter("raw", new double[][]{
                {3.0, 0.0},
                {4.0, 5.0}
        });
        final NormalizedMatrixParameter matrix = new NormalizedMatrixParameter(
                "normalized", raw, NormalizedMatrixParameter.Normalization.COLUMNS);

        assertEquals(5.0, matrix.getColumnNorm(0), 1e-15);
        assertEquals(5.0, matrix.getColumnNorm(1), 1e-15);
        assertEquals(0.6, matrix.getParameterValue(0, 0), 1e-15);
        assertEquals(0.8, matrix.getParameterValue(1, 0), 1e-15);
        assertEquals(0.0, matrix.getParameterValue(0, 1), 1e-15);
        assertEquals(1.0, matrix.getParameterValue(1, 1), 1e-15);
    }

    public void testGlobalPullBackGradientMatchesFiniteDifference() {
        assertPullBackMatchesFiniteDifference(NormalizedMatrixParameter.Normalization.GLOBAL);
    }

    public void testColumnPullBackGradientMatchesFiniteDifference() {
        assertPullBackMatchesFiniteDifference(NormalizedMatrixParameter.Normalization.COLUMNS);
    }

    public void testCanonicalProjectorRoutesGradientToRawMatrix() {
        final MatrixParameter raw = makeMatrixParameter("raw", new double[][]{
                {2.0, 1.0},
                {0.5, 3.0}
        });
        final NormalizedMatrixParameter rotation = new NormalizedMatrixParameter(
                "normalized", raw, NormalizedMatrixParameter.Normalization.GLOBAL);
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
        final double[] expectedRotation = rotation.pullBackGradientForParameter(raw, gradientR, 2);

        final double[] rawGradient = CanonicalSelectionGradientProjector.assembleBlockGradientResultFlat(
                2, raw, block, nativeBlockGradient, gradientR);
        for (int i = 0; i < expectedRotation.length; i++) {
            assertEquals(expectedRotation[i], rawGradient[i], 1e-15);
        }

        final double[] compoundGradient = CanonicalSelectionGradientProjector.assembleBlockGradientResultFlat(
                2, block.getParameter(), block, nativeBlockGradient, gradientR);
        assertEquals(nativeBlockGradient.length + expectedRotation.length, compoundGradient.length);
        for (int i = 0; i < nativeBlockGradient.length; i++) {
            assertEquals(nativeBlockGradient[i], compoundGradient[i], 1e-15);
        }
        for (int i = 0; i < expectedRotation.length; i++) {
            assertEquals(expectedRotation[i], compoundGradient[nativeBlockGradient.length + i], 1e-15);
        }
    }

    public void testGlobalScalingPreservesDenseDriftWithFixedBlock() {
        final MatrixParameter raw = makeMatrixParameter("raw", new double[][]{
                {1.0, 0.3, -0.4},
                {0.2, 1.4, 0.5},
                {-0.1, 0.7, 1.2}
        });
        final MatrixParameter scaled = makeMatrixParameter("scaled", new double[][]{
                {7.0, 2.1, -2.8},
                {1.4, 9.8, 3.5},
                {-0.7, 4.9, 8.4}
        });
        final NormalizedMatrixParameter left = new NormalizedMatrixParameter(
                "left", raw, NormalizedMatrixParameter.Normalization.GLOBAL);
        final NormalizedMatrixParameter right = new NormalizedMatrixParameter(
                "right", scaled, NormalizedMatrixParameter.Normalization.GLOBAL);
        final double[][] block = new double[][]{
                {1.0, 0.0, 0.0},
                {0.0, 0.8, -0.2},
                {0.0, 0.2, 0.8}
        };

        assertMatrixEquals(drift(left, block), drift(right, block), 1e-13);
    }

    public void testRejectsNearZeroNorms() {
        final MatrixParameter raw = makeMatrixParameter("raw", new double[][]{
                {1.0, 0.0},
                {0.0, 0.0}
        });
        try {
            new NormalizedMatrixParameter(
                    "normalized", raw, NormalizedMatrixParameter.Normalization.COLUMNS, 1e-10)
                    .getParameterValue(0, 1);
            fail("Expected IllegalArgumentException for near-zero raw column");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage().contains("Column 1"));
        }

        final MatrixParameter zero = makeMatrixParameter("zero", new double[][]{
                {0.0, 0.0},
                {0.0, 0.0}
        });
        try {
            new NormalizedMatrixParameter(
                    "normalized", zero, NormalizedMatrixParameter.Normalization.GLOBAL, 1e-10)
                    .getParameterValue(0, 0);
            fail("Expected IllegalArgumentException for near-zero raw matrix");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage().contains("Frobenius norm"));
        }
    }

    private static void assertPullBackMatchesFiniteDifference(
            final NormalizedMatrixParameter.Normalization normalization) {
        final MatrixParameter raw = makeMatrixParameter("raw", new double[][]{
                {1.3, -0.4, 0.7},
                {0.2, 1.1, -0.6},
                {-0.5, 0.8, 1.5}
        });
        final NormalizedMatrixParameter matrix =
                new NormalizedMatrixParameter("normalized", raw, normalization);
        final double[] gradientR = new double[]{
                0.7, -1.1, 0.4,
                0.3, 0.9, -0.8,
                -0.2, 1.4, 0.6
        };

        final double[] analytic = matrix.pullBackGradientForParameter(raw, gradientR, 3);
        for (int i = 0; i < analytic.length; i++) {
            assertEquals(
                    "finite difference gradient for raw[" + i + "]",
                    finiteDifference(matrix, raw, i, gradientR),
                    analytic[i],
                    1e-6);
        }
    }

    private static MatrixParameter makeMatrixParameter(final String name,
                                                       final double[][] values) {
        final int rows = values.length;
        final int cols = values[0].length;
        final Parameter[] columns = new Parameter[cols];
        for (int col = 0; col < cols; col++) {
            final double[] column = new double[rows];
            for (int row = 0; row < rows; row++) {
                column[row] = values[row][col];
            }
            columns[col] = new Parameter.Default(name + ".col" + col, column);
        }
        return new MatrixParameter(name, columns);
    }

    private static double finiteDifference(final NormalizedMatrixParameter matrix,
                                           final MatrixParameter raw,
                                           final int nativeIndex,
                                           final double[] gradientR) {
        final int dim = matrix.getRowDimension();
        final int row = nativeIndex % dim;
        final int col = nativeIndex / dim;
        final double h = 1e-7;
        final double original = raw.getParameterValue(row, col);
        raw.setParameterValue(row, col, original + h);
        matrix.fireParameterChangedEvent();
        final double plus = objective(matrix, gradientR);
        raw.setParameterValue(row, col, original - h);
        matrix.fireParameterChangedEvent();
        final double minus = objective(matrix, gradientR);
        raw.setParameterValue(row, col, original);
        matrix.fireParameterChangedEvent();
        return (plus - minus) / (2.0 * h);
    }

    private static double objective(final NormalizedMatrixParameter matrix,
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

    private static double frobeniusNorm(final NormalizedMatrixParameter matrix) {
        double sum = 0.0;
        final int dim = matrix.getRowDimension();
        for (int row = 0; row < dim; row++) {
            for (int col = 0; col < dim; col++) {
                final double value = matrix.getParameterValue(row, col);
                sum += value * value;
            }
        }
        return Math.sqrt(sum);
    }

    private static double[][] drift(final NormalizedMatrixParameter rotation,
                                    final double[][] block) {
        final int dim = rotation.getRowDimension();
        final double[][] r = rotation.getParameterAsMatrix();
        final double[][] inverse = invert3(r);
        return multiply(multiply(r, block), inverse);
    }

    private static double[][] multiply(final double[][] left, final double[][] right) {
        final int dim = left.length;
        final double[][] out = new double[dim][dim];
        for (int row = 0; row < dim; row++) {
            for (int col = 0; col < dim; col++) {
                double sum = 0.0;
                for (int k = 0; k < dim; k++) {
                    sum += left[row][k] * right[k][col];
                }
                out[row][col] = sum;
            }
        }
        return out;
    }

    private static double[][] invert3(final double[][] m) {
        final double determinant =
                m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
                        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
                        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
        final double[][] out = new double[3][3];
        out[0][0] = (m[1][1] * m[2][2] - m[1][2] * m[2][1]) / determinant;
        out[0][1] = (m[0][2] * m[2][1] - m[0][1] * m[2][2]) / determinant;
        out[0][2] = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) / determinant;
        out[1][0] = (m[1][2] * m[2][0] - m[1][0] * m[2][2]) / determinant;
        out[1][1] = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) / determinant;
        out[1][2] = (m[0][2] * m[1][0] - m[0][0] * m[1][2]) / determinant;
        out[2][0] = (m[1][0] * m[2][1] - m[1][1] * m[2][0]) / determinant;
        out[2][1] = (m[0][1] * m[2][0] - m[0][0] * m[2][1]) / determinant;
        out[2][2] = (m[0][0] * m[1][1] - m[0][1] * m[1][0]) / determinant;
        return out;
    }

    private static void assertMatrixEquals(final double[][] expected,
                                           final double[][] actual,
                                           final double tolerance) {
        for (int row = 0; row < expected.length; row++) {
            for (int col = 0; col < expected.length; col++) {
                assertEquals(expected[row][col], actual[row][col], tolerance);
            }
        }
    }
}
