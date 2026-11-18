package test.dr.evomodel.treedatalikelihood.discrete;

import dr.evomodel.treedatalikelihood.discrete.FisherWhiteningBuilder;
import junit.framework.TestCase;

public final class FisherWhiteningBuilderTest extends TestCase {

    public void testRankReducedWhiteningUsesPriorInNullSpace() {
        final double[][] fisher = new double[][]{
                {100.0, 0.0, 0.0},
                {0.0, 4.0, 0.0},
                {0.0, 0.0, 0.0}
        };
        final double[][] prior = new double[][]{
                {1.0, 0.0, 0.0},
                {0.0, 1.0, 0.0},
                {0.0, 0.0, 9.0}
        };
        final FisherWhiteningBuilder.Result result = FisherWhiteningBuilder.build(
                fisher, prior, 10.0, 1.0e-6, Double.POSITIVE_INFINITY);

        assertEquals(2, result.rank);
        assertEquals(1, result.nullity);
        assertEquals(4.0, result.minimumPositiveEigenvalue, 1.0e-12);
        assertEquals(100.0, result.maximumEigenvalue, 1.0e-12);

        final double[][] covariance = multiplyByTranspose(result.map);
        // Median-positive/fraction is 10 here, so the weaker identifiable
        // direction is floored from precision 4 to precision 10.
        assertEquals(1.0 / 100.0, covariance[0][0], 1.0e-10);
        assertEquals(1.0 / 10.0, covariance[1][1], 1.0e-10);
        assertEquals(1.0 / 9.0, covariance[2][2], 1.0e-10);
        assertEquals(0.0, covariance[0][1], 1.0e-10);
        assertEquals(0.0, covariance[0][2], 1.0e-10);
        assertEquals(0.0, covariance[1][2], 1.0e-10);
    }

    public void testTrustRadiusCapsWhiteningColumns() {
        final double[][] fisher = new double[][]{
                {1.0e-4, 0.0},
                {0.0, 1.0}
        };
        final FisherWhiteningBuilder.Result result = FisherWhiteningBuilder.build(
                fisher, null, 1.0e12, 1.0, 0.25);
        for (int column = 0; column < 2; ++column) {
            double squaredNorm = 0.0;
            for (int row = 0; row < 2; ++row) {
                squaredNorm += result.map[row][column] * result.map[row][column];
            }
            assertTrue(Math.sqrt(squaredNorm) <= 0.25 + 1.0e-12);
        }
    }

    private static double[][] multiplyByTranspose(final double[][] matrix) {
        final double[][] result = new double[matrix.length][matrix.length];
        for (int i = 0; i < matrix.length; ++i) {
            for (int j = 0; j < matrix.length; ++j) {
                for (int k = 0; k < matrix.length; ++k) {
                    result[i][j] += matrix[i][k] * matrix[j][k];
                }
            }
        }
        return result;
    }
}
