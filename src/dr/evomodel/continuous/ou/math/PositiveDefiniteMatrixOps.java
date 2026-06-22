package dr.evomodel.continuous.ou.math;

import org.ejml.data.DenseMatrix64F;

/**
 * Allocation-free positive-definite matrix operations for OU branch kernels.
 */
public final class PositiveDefiniteMatrixOps {

    private PositiveDefiniteMatrixOps() {
    }

    public static double copyAndInvertPositiveDefiniteFlat(
            final DenseMatrix64F source,
            final double[] matrixOut,
            final double[] inverseOut,
            final double[] choleskyScratch,
            final double[] lowerInverseScratch,
            final PositiveDefiniteMatrixInversionOptions options) {
        final int d = source.numRows;
        final double[] sourceData = source.data;
        double maxAbsDiagonal = 0.0;
        for (int i = 0; i < d; ++i) {
            final int iOff = i * d;
            for (int j = 0; j < d; ++j) {
                final double value = 0.5 * (sourceData[iOff + j] + sourceData[j * d + i]);
                matrixOut[iOff + j] = value;
                choleskyScratch[iOff + j] = value;
            }
            maxAbsDiagonal = Math.max(maxAbsDiagonal, Math.abs(matrixOut[iOff + i]));
        }

        return invertSymmetrizedPositiveDefiniteFlat(
                inverseOut,
                choleskyScratch,
                lowerInverseScratch,
                d,
                maxAbsDiagonal,
                options);
    }

    public static double copyAndInvertPositiveDefiniteFlat(
            final double[] source,
            final int dimension,
            final double[] matrixOut,
            final double[] inverseOut,
            final double[] choleskyScratch,
            final double[] lowerInverseScratch,
            final PositiveDefiniteMatrixInversionOptions options) {
        final int d = dimension;
        double maxAbsDiagonal = 0.0;
        for (int i = 0; i < d; ++i) {
            final int iOff = i * d;
            for (int j = 0; j < d; ++j) {
                final double value = 0.5 * (source[iOff + j] + source[j * d + i]);
                matrixOut[iOff + j] = value;
                choleskyScratch[iOff + j] = value;
            }
            maxAbsDiagonal = Math.max(maxAbsDiagonal, Math.abs(matrixOut[iOff + i]));
        }

        return invertSymmetrizedPositiveDefiniteFlat(
                inverseOut,
                choleskyScratch,
                lowerInverseScratch,
                d,
                maxAbsDiagonal,
                options);
    }

    private static double invertSymmetrizedPositiveDefiniteFlat(
            final double[] inverseOut,
            final double[] choleskyScratch,
            final double[] lowerInverseScratch,
            final int dimension,
            final double maxAbsDiagonal,
            final PositiveDefiniteMatrixInversionOptions options) {
        final double pivotFloor = options.jitterBase(maxAbsDiagonal);
        final PivotStats pivotStats = factorCholeskyFlatInPlace(
                choleskyScratch, dimension, pivotFloor);
        if (options.isPivotFloorDebugEnabled()) {
            options.reportPivotFloor(
                    pivotStats.clippedPivotCount,
                    dimension,
                    pivotStats.minPivotBeforeFloor,
                    pivotFloor,
                    maxAbsDiagonal);
        }
        invertFlatCholesky(choleskyScratch, lowerInverseScratch, dimension, pivotFloor);

        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                double sum = 0.0;
                for (int k = 0; k < dimension; ++k) {
                    sum += lowerInverseScratch[k * dimension + i]
                            * lowerInverseScratch[k * dimension + j];
                }
                inverseOut[i * dimension + j] = sum;
            }
        }
        return 2.0 * pivotStats.logDet;
    }

    private static PivotStats factorCholeskyFlatInPlace(final double[] cholesky,
                                                        final int dimension,
                                                        final double pivotFloor) {
        int clippedPivotCount = 0;
        double minPivotBeforeFloor = Double.POSITIVE_INFINITY;
        double logDet = 0.0;
        for (int i = 0; i < dimension; ++i) {
            final int iOff = i * dimension;
            for (int j = 0; j <= i; ++j) {
                double sum = cholesky[iOff + j];
                for (int k = 0; k < j; ++k) {
                    sum -= cholesky[iOff + k] * cholesky[j * dimension + k];
                }
                if (i == j) {
                    minPivotBeforeFloor = Math.min(minPivotBeforeFloor, sum);
                    if (sum < pivotFloor) {
                        clippedPivotCount++;
                        sum = pivotFloor;
                    }
                    cholesky[iOff + i] = Math.sqrt(sum);
                    logDet += Math.log(cholesky[iOff + i]);
                } else {
                    final double denominator =
                            Math.max(cholesky[j * dimension + j], Math.sqrt(pivotFloor));
                    cholesky[iOff + j] = sum / denominator;
                }
            }
        }
        return new PivotStats(clippedPivotCount, minPivotBeforeFloor, logDet);
    }

    private static void invertFlatCholesky(final double[] cholesky,
                                           final double[] lowerInverse,
                                           final int dimension,
                                           final double pivotFloor) {
        for (int col = 0; col < dimension; ++col) {
            for (int row = 0; row < dimension; ++row) {
                double sum = (row == col) ? 1.0 : 0.0;
                final int rowOff = row * dimension;
                for (int k = 0; k < row; ++k) {
                    sum -= cholesky[rowOff + k] * lowerInverse[k * dimension + col];
                }
                final double denominator =
                        Math.max(cholesky[rowOff + row], Math.sqrt(pivotFloor));
                lowerInverse[rowOff + col] = sum / denominator;
            }
        }
    }

    private static final class PivotStats {
        private final int clippedPivotCount;
        private final double minPivotBeforeFloor;
        private final double logDet;

        private PivotStats(final int clippedPivotCount,
                           final double minPivotBeforeFloor,
                           final double logDet) {
            this.clippedPivotCount = clippedPivotCount;
            this.minPivotBeforeFloor = minPivotBeforeFloor;
            this.logDet = logDet;
        }
    }
}
