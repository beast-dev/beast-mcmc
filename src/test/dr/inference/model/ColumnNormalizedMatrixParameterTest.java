package test.dr.inference.model;

import dr.evomodel.treedatalikelihood.continuous.canonical.gradient.CanonicalSelectionGradientProjector;
import dr.inference.model.BlockDiagonalPolarStableMatrixParameter;
import dr.inference.model.ColumnNormalizedMatrixParameter;
import dr.inference.model.MatrixParameter;
import dr.inference.model.NativeMatrixGradientProvider;
import dr.inference.model.Parameter;
import test.dr.math.MathTestCase;

public class ColumnNormalizedMatrixParameterTest extends MathTestCase {

    public void testColumnsHaveUnitNormAndExposeRawMatrixAsNativeParameter() {
        final MatrixParameter raw = makeMatrixParameter("raw", new double[][]{
                {3.0, 0.0},
                {4.0, 5.0}
        });
        final ColumnNormalizedMatrixParameter matrix =
                new ColumnNormalizedMatrixParameter("normalized", raw);

        assertTrue(matrix instanceof NativeMatrixGradientProvider);
        assertTrue(matrix.supportsNativeParameter(raw));
        assertTrue(matrix.getNativeParameter() == raw);
        assertEquals(4, matrix.getPullBackGradientDimension(raw));
        assertEquals(5.0, matrix.getColumnNorm(0), 1e-15);
        assertEquals(5.0, matrix.getColumnNorm(1), 1e-15);

        assertEquals(0.6, matrix.getParameterValue(0, 0), 1e-15);
        assertEquals(0.8, matrix.getParameterValue(1, 0), 1e-15);
        assertEquals(0.0, matrix.getParameterValue(0, 1), 1e-15);
        assertEquals(1.0, matrix.getParameterValue(1, 1), 1e-15);
    }

    public void testPullBackGradientMatchesFiniteDifference() {
        final MatrixParameter raw = makeMatrixParameter("raw", new double[][]{
                {1.3, -0.4, 0.7},
                {0.2, 1.1, -0.6},
                {-0.5, 0.8, 1.5}
        });
        final ColumnNormalizedMatrixParameter matrix =
                new ColumnNormalizedMatrixParameter("normalized", raw);
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

    public void testCanonicalProjectorRoutesGradientToRawMatrix() {
        final MatrixParameter raw = makeMatrixParameter("raw", new double[][]{
                {2.0, 1.0},
                {0.5, 3.0}
        });
        final ColumnNormalizedMatrixParameter rotation =
                new ColumnNormalizedMatrixParameter("normalized", raw);
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

    public void testRejectsNearZeroColumnNorm() {
        final MatrixParameter raw = makeMatrixParameter("raw", new double[][]{
                {1.0, 0.0},
                {0.0, 0.0}
        });
        final ColumnNormalizedMatrixParameter matrix =
                new ColumnNormalizedMatrixParameter("normalized", raw, 1e-10);

        try {
            matrix.getParameterValue(0, 1);
            fail("Expected IllegalArgumentException for near-zero raw column");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage().contains("Column 1"));
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

    private static double finiteDifference(final ColumnNormalizedMatrixParameter matrix,
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

    private static double objective(final ColumnNormalizedMatrixParameter matrix,
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
