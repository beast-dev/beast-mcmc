package dr.evomodel.continuous.ou.orthogonalblockdiagonal;

import org.ejml.data.DenseMatrix64F;

/**
 * Positive-definite inversion helper used by orthogonal block canonical OU code.
 */
final class OrthogonalBlockPositiveDefiniteInverter {

    private static final double SPD_JITTER_RELATIVE = 1.0e-14;
    private static final double SPD_JITTER_ABSOLUTE = 1.0e-14;
    private static final String DEBUG_PD_PIVOT_FLOOR_PROPERTY =
            "beast.debug.pdPivotFloor";

    private OrthogonalBlockPositiveDefiniteInverter() { }

    static double copyAndInvert(final DenseMatrix64F source,
                                final double[][] matrixOut,
                                final double[][] inverseOut,
                                final double[] choleskyScratch,
                                final double[] lowerInverseScratch) {
        final int d = source.numRows;
        final double[] sourceData = source.data;
        double maxAbsDiagonal = 0.0;
        for (int i = 0; i < d; ++i) {
            final double[] matrixRow = matrixOut[i];
            final int iOff = i * d;
            for (int j = 0; j < d; ++j) {
                final double value = 0.5 * (sourceData[iOff + j] + sourceData[j * d + i]);
                matrixRow[j] = value;
                choleskyScratch[iOff + j] = value;
            }
            maxAbsDiagonal = Math.max(maxAbsDiagonal, Math.abs(matrixRow[i]));
        }

        final double pivotFloor = pivotFloor(maxAbsDiagonal);
        final PivotStats pivotStats = factorCholeskyInPlace(choleskyScratch, d, pivotFloor);
        maybeReportPivotFloor(pivotStats, d, pivotFloor, maxAbsDiagonal);
        invertFromCholesky(choleskyScratch, lowerInverseScratch, d, pivotFloor);

        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                double sum = 0.0;
                for (int k = 0; k < d; ++k) {
                    sum += lowerInverseScratch[k * d + i] * lowerInverseScratch[k * d + j];
                }
                inverseOut[i][j] = sum;
            }
        }
        return 2.0 * pivotStats.logDet;
    }

    static double copyAndInvertFlat(final DenseMatrix64F source,
                                    final double[] matrixOut,
                                    final double[] inverseOut,
                                    final double[] choleskyScratch,
                                    final double[] lowerInverseScratch) {
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

        final double pivotFloor = pivotFloor(maxAbsDiagonal);
        final PivotStats pivotStats = factorCholeskyInPlace(choleskyScratch, d, pivotFloor);
        maybeReportPivotFloor(pivotStats, d, pivotFloor, maxAbsDiagonal);
        invertFromCholesky(choleskyScratch, lowerInverseScratch, d, pivotFloor);

        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                double sum = 0.0;
                for (int k = 0; k < d; ++k) {
                    sum += lowerInverseScratch[k * d + i] * lowerInverseScratch[k * d + j];
                }
                inverseOut[i * d + j] = sum;
            }
        }
        return 2.0 * pivotStats.logDet;
    }

    private static double pivotFloor(final double maxAbsDiagonal) {
        return Math.max(
                SPD_JITTER_ABSOLUTE,
                SPD_JITTER_RELATIVE * Math.max(1.0, maxAbsDiagonal));
    }

    private static PivotStats factorCholeskyInPlace(final double[] cholesky,
                                                    final int d,
                                                    final double pivotFloor) {
        int clippedPivotCount = 0;
        double minPivotBeforeFloor = Double.POSITIVE_INFINITY;
        double logDet = 0.0;
        for (int i = 0; i < d; ++i) {
            final int iOff = i * d;
            for (int j = 0; j <= i; ++j) {
                double sum = cholesky[iOff + j];
                for (int k = 0; k < j; ++k) {
                    sum -= cholesky[iOff + k] * cholesky[j * d + k];
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
                    final double denominator = Math.max(cholesky[j * d + j], Math.sqrt(pivotFloor));
                    cholesky[iOff + j] = sum / denominator;
                }
            }
        }
        return new PivotStats(clippedPivotCount, minPivotBeforeFloor, logDet);
    }

    private static void invertFromCholesky(final double[] cholesky,
                                           final double[] lowerInverse,
                                           final int d,
                                           final double pivotFloor) {
        for (int col = 0; col < d; ++col) {
            for (int row = 0; row < d; ++row) {
                double sum = (row == col) ? 1.0 : 0.0;
                final int rowOff = row * d;
                for (int k = 0; k < row; ++k) {
                    sum -= cholesky[rowOff + k] * lowerInverse[k * d + col];
                }
                final double denominator = Math.max(cholesky[rowOff + row], Math.sqrt(pivotFloor));
                lowerInverse[rowOff + col] = sum / denominator;
            }
        }
    }

    private static void maybeReportPivotFloor(final PivotStats pivotStats,
                                              final int dimension,
                                              final double pivotFloor,
                                              final double maxAbsDiagonal) {
        if (Boolean.getBoolean(DEBUG_PD_PIVOT_FLOOR_PROPERTY) && pivotStats.clippedPivotCount > 0) {
            System.err.println(
                    "pdPivotFloorDebug clipped=" + pivotStats.clippedPivotCount
                            + " dim=" + dimension
                            + " minPivotBeforeFloor=" + pivotStats.minPivotBeforeFloor
                            + " pivotFloor=" + pivotFloor
                            + " maxAbsDiagonal=" + maxAbsDiagonal);
        }
    }

    private static final class PivotStats {
        final int clippedPivotCount;
        final double minPivotBeforeFloor;
        final double logDet;

        private PivotStats(final int clippedPivotCount,
                           final double minPivotBeforeFloor,
                           final double logDet) {
            this.clippedPivotCount = clippedPivotCount;
            this.minPivotBeforeFloor = minPivotBeforeFloor;
            this.logDet = logDet;
        }
    }
}
