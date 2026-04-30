package dr.evomodel.treedatalikelihood.continuous.gaussian.message;

import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;

/**
 * Shared numerical kernels for canonical Gaussian machinery.
 *
 * <p>The methods here deliberately operate on caller-owned buffers. That keeps
 * hot likelihood and gradient paths allocation-free while giving the canonical
 * code one place for common stability policy.</p>
 */
public final class CanonicalNumerics {

    private static final double SYMMETRIC_JITTER_RELATIVE = 1.0e-12;
    private static final double SYMMETRIC_JITTER_ABSOLUTE = 1.0e-12;
    private static final int ROBUST_INVERSION_ATTEMPTS = 12;

    private CanonicalNumerics() { }

    public interface DenseSpdFailureDump {
        void appendJson(StringBuilder sb,
                        String context,
                        DenseMatrix64F originalSource,
                        DenseMatrix64F symmetrizedSource,
                        double jitterBase);
    }

    /**
     * Robustly invert a symmetric positive-definite matrix.
     *
     * <p>The source is symmetrized into {@code symmetricScratch}; retries add
     * diagonal jitter to {@code choleskyScratch}. The returned value is the
     * log-determinant of the original matrix after any applied jitter.</p>
     */
    public static double invertSymmetricPositiveDefinite(
            final double[][] matrix,
            final int dimension,
            final double[][] inverseOut,
            final double[][] symmetricScratch,
            final double[][] choleskyScratch,
            final double[][] lowerInverseScratch) {
        checkSquare(matrix, dimension, "matrix");
        checkSquare(inverseOut, dimension, "inverseOut");
        checkSquare(symmetricScratch, dimension, "symmetricScratch");
        checkSquare(choleskyScratch, dimension, "choleskyScratch");
        checkSquare(lowerInverseScratch, dimension, "lowerInverseScratch");

        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                symmetricScratch[i][j] = 0.5 * (matrix[i][j] + matrix[j][i]);
            }
        }

        final double jitterBase = Math.max(
                SYMMETRIC_JITTER_ABSOLUTE,
                SYMMETRIC_JITTER_RELATIVE * Math.max(1.0, maxAbsDiagonal(symmetricScratch, dimension)));
        final double lowerBound = gershgorinLowerBound(symmetricScratch, dimension);

        double jitter = 0.0;
        for (int attempt = 0; attempt < ROBUST_INVERSION_ATTEMPTS; ++attempt) {
            copySquare(symmetricScratch, choleskyScratch, dimension);
            if (jitter > 0.0) {
                addDiagonal(choleskyScratch, dimension, jitter);
            }

            final double logDet = invertCholeskyInPlace(
                    choleskyScratch, lowerInverseScratch, inverseOut, dimension);
            if (Double.isFinite(logDet) && isFinite(inverseOut, dimension)) {
                symmetrizeInPlace(inverseOut, dimension);
                return logDet;
            }

            if (jitter == 0.0) {
                jitter = lowerBound > 0.0 ? jitterBase : (-lowerBound + jitterBase);
            } else {
                jitter *= 10.0;
            }
        }

        throw new IllegalStateException(
                "Failed to invert symmetric positive definite matrix stably"
                        + "; dim=" + dimension
                        + "; finite=" + isFinite(symmetricScratch, dimension)
                        + "; gershgorinLowerBound=" + gershgorinLowerBound(symmetricScratch, dimension)
                        + "; minDiag=" + minDiagonal(symmetricScratch, dimension)
                        + "; maxAbsDiag=" + maxAbsDiagonal(symmetricScratch, dimension));
    }

    public static void safeInvertSymmetricPositiveDefinite(
            final DenseMatrix64F source,
            final DenseMatrix64F inverseOut,
            final DenseMatrix64F symmetricScratch,
            final DenseMatrix64F workingScratch,
            final double[][] symmetric2DScratch,
            final double[][] adjusted2DScratch,
            final double[][] choleskyScratch,
            final double[][] lowerInverseScratch,
            final CanonicalNumericsOptions options,
            final String context,
            final DenseSpdFailureDump failureDump) {
        symmetrizeCopy(source, symmetricScratch);
        final double jitterBase = options.jitterBase(maxAbsDiagonal(symmetricScratch));

        if (options.isForceStrictSpdInversion()) {
            if (invertSymmetricPositiveDefiniteStrictFallback(
                    symmetricScratch,
                    inverseOut,
                    jitterBase,
                    options.getStrictInversionAttempts(),
                    symmetric2DScratch,
                    adjusted2DScratch,
                    choleskyScratch,
                    lowerInverseScratch)) {
                return;
            }
            emitSpdFailureDebugDump(options, context, source, symmetricScratch, jitterBase, failureDump);
            throw new IllegalStateException(
                    "Failed strict SPD inversion with fallback; context=" + context);
        }

        double jitter = 0.0;
        for (int attempt = 0; attempt < options.getRobustInversionAttempts(); ++attempt) {
            workingScratch.set(symmetricScratch);
            if (jitter > 0.0) {
                addDiagonalJitter(workingScratch, jitter);
            }
            if (CommonOps.invert(workingScratch, inverseOut) && isFinite(inverseOut)) {
                symmetrizeInPlace(inverseOut);
                return;
            }
            jitter = jitter == 0.0 ? jitterBase : 10.0 * jitter;
        }
        if (invertSymmetricPositiveDefiniteStrictFallback(
                symmetricScratch,
                inverseOut,
                jitterBase,
                options.getStrictInversionAttempts(),
                symmetric2DScratch,
                adjusted2DScratch,
                choleskyScratch,
                lowerInverseScratch)) {
            return;
        }
        emitSpdFailureDebugDump(options, context, source, symmetricScratch, jitterBase, failureDump);
        throw new IllegalStateException(
                "Failed to invert symmetric positive definite matrix stably; context=" + context);
    }

    public static void safeInvertWithJitter(final DenseMatrix64F source,
                                            final DenseMatrix64F inverseOut,
                                            final DenseMatrix64F workingScratch,
                                            final CanonicalNumericsOptions options) {
        final double jitterBase = options.jitterBase(maxAbsDiagonal(source));
        double jitter = jitterBase;
        for (int attempt = 0; attempt < options.getRobustInversionAttempts(); ++attempt) {
            workingScratch.set(source);
            addDiagonalJitter(workingScratch, jitter);
            if (CommonOps.invert(workingScratch, inverseOut) && isFinite(inverseOut)) {
                return;
            }
            jitter *= 10.0;
        }
        throw new IllegalStateException("Failed to invert matrix stably");
    }

    public static double copyAndInvertPositiveDefinite(
            final DenseMatrix64F source,
            final double[][] matrixOut,
            final double[][] inverseOut,
            final double[] choleskyScratch,
            final double[] lowerInverseScratch,
            final CanonicalNumericsOptions options) {
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

        final double pivotFloor = options.jitterBase(maxAbsDiagonal);
        final PivotStats pivotStats = factorCholeskyFlatInPlace(choleskyScratch, d, pivotFloor);
        maybeReportPivotFloor(options, pivotStats, d, pivotFloor, maxAbsDiagonal);
        invertFlatCholesky(choleskyScratch, lowerInverseScratch, d, pivotFloor);

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

    public static double copyAndInvertPositiveDefiniteFlat(
            final DenseMatrix64F source,
            final double[] matrixOut,
            final double[] inverseOut,
            final double[] choleskyScratch,
            final double[] lowerInverseScratch,
            final CanonicalNumericsOptions options) {
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

        final double pivotFloor = options.jitterBase(maxAbsDiagonal);
        final PivotStats pivotStats = factorCholeskyFlatInPlace(choleskyScratch, d, pivotFloor);
        maybeReportPivotFloor(options, pivotStats, d, pivotFloor, maxAbsDiagonal);
        invertFlatCholesky(choleskyScratch, lowerInverseScratch, d, pivotFloor);

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

    private static double invertCholeskyInPlace(final double[][] cholesky,
                                                final double[][] lowerInverse,
                                                final double[][] inverseOut,
                                                final int dimension) {
        double logDet = 0.0;
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j <= i; ++j) {
                double sum = cholesky[i][j];
                for (int k = 0; k < j; ++k) {
                    sum -= cholesky[i][k] * cholesky[j][k];
                }
                if (i == j) {
                    if (!(sum > 0.0) || !Double.isFinite(sum)) {
                        return Double.NaN;
                    }
                    final double diag = Math.sqrt(sum);
                    cholesky[i][j] = diag;
                    logDet += Math.log(diag);
                } else {
                    cholesky[i][j] = sum / cholesky[j][j];
                }
            }
            for (int j = i + 1; j < dimension; ++j) {
                cholesky[i][j] = 0.0;
            }
        }
        logDet *= 2.0;

        for (int column = 0; column < dimension; ++column) {
            for (int row = 0; row < dimension; ++row) {
                double sum = row == column ? 1.0 : 0.0;
                for (int k = 0; k < row; ++k) {
                    sum -= cholesky[row][k] * lowerInverse[k][column];
                }
                lowerInverse[row][column] = sum / cholesky[row][row];
            }
        }

        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                double sum = 0.0;
                for (int k = 0; k < dimension; ++k) {
                    sum += lowerInverse[k][i] * lowerInverse[k][j];
                }
                inverseOut[i][j] = sum;
            }
        }
        return logDet;
    }

    private static boolean invertSymmetricPositiveDefiniteStrictFallback(
            final DenseMatrix64F symmetricSource,
            final DenseMatrix64F inverseOut,
            final double jitterBase,
            final int attempts,
            final double[][] symmetric,
            final double[][] adjusted,
            final double[][] cholesky,
            final double[][] lowerInverse) {
        final int d = symmetricSource.numRows;
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                symmetric[i][j] = symmetricSource.unsafe_get(i, j);
            }
        }

        double jitter = 0.0;
        for (int attempt = 0; attempt < attempts; ++attempt) {
            copySquare(symmetric, adjusted, d);
            if (jitter > 0.0) {
                addDiagonal(adjusted, d, jitter);
            }
            if (invertSymmetricPositiveDefiniteStrict(adjusted, cholesky, lowerInverse, inverseOut, d)) {
                symmetrizeInPlace(inverseOut);
                return true;
            }
            jitter = jitter == 0.0 ? jitterBase : 10.0 * jitter;
        }
        return false;
    }

    private static boolean invertSymmetricPositiveDefiniteStrict(final double[][] matrix,
                                                                 final double[][] cholesky,
                                                                 final double[][] lowerInverse,
                                                                 final DenseMatrix64F inverseOut,
                                                                 final int dimensionUsed) {
        final int d = dimensionUsed;
        copySquare(matrix, cholesky, d);

        for (int i = 0; i < d; ++i) {
            for (int j = 0; j <= i; ++j) {
                double sum = cholesky[i][j];
                for (int k = 0; k < j; ++k) {
                    sum -= cholesky[i][k] * cholesky[j][k];
                }
                if (i == j) {
                    if (!(sum > 0.0) || !Double.isFinite(sum)) {
                        return false;
                    }
                    cholesky[i][j] = Math.sqrt(sum);
                } else {
                    cholesky[i][j] = sum / cholesky[j][j];
                }
            }
            for (int j = i + 1; j < d; ++j) {
                cholesky[i][j] = 0.0;
            }
        }

        for (int column = 0; column < d; ++column) {
            for (int row = 0; row < d; ++row) {
                double sum = row == column ? 1.0 : 0.0;
                for (int k = 0; k < row; ++k) {
                    sum -= cholesky[row][k] * lowerInverse[k][column];
                }
                lowerInverse[row][column] = sum / cholesky[row][row];
            }
        }

        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                double sum = 0.0;
                for (int k = 0; k < d; ++k) {
                    sum += lowerInverse[k][i] * lowerInverse[k][j];
                }
                inverseOut.unsafe_set(i, j, sum);
            }
        }
        return isFinite(inverseOut);
    }

    private static PivotStats factorCholeskyFlatInPlace(final double[] cholesky,
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

    private static void invertFlatCholesky(final double[] cholesky,
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

    private static void maybeReportPivotFloor(final CanonicalNumericsOptions options,
                                              final PivotStats pivotStats,
                                              final int dimension,
                                              final double pivotFloor,
                                              final double maxAbsDiagonal) {
        if (options.isPivotFloorDebugEnabled() && pivotStats.clippedPivotCount > 0) {
            System.err.println(
                    "pdPivotFloorDebug clipped=" + pivotStats.clippedPivotCount
                            + " dim=" + dimension
                            + " minPivotBeforeFloor=" + pivotStats.minPivotBeforeFloor
                            + " pivotFloor=" + pivotFloor
                            + " maxAbsDiagonal=" + maxAbsDiagonal);
        }
    }

    private static void emitSpdFailureDebugDump(final CanonicalNumericsOptions options,
                                                final String context,
                                                final DenseMatrix64F originalSource,
                                                final DenseMatrix64F symmetrizedSource,
                                                final double jitterBase,
                                                final DenseSpdFailureDump failureDump) {
        if (!options.isSpdFailureDumpEnabled()) {
            return;
        }
        final String path = "/tmp/ou_spd_failure_" + System.nanoTime() + ".json";
        try {
            final StringBuilder sb = new StringBuilder(8192);
            sb.append("{\n");
            sb.append("\"context\":\"").append(context).append("\",\n");
            sb.append("\"dimension\":").append(symmetrizedSource.numRows).append(",\n");
            sb.append("\"jitterBase\":").append(jitterBase).append(",\n");
            sb.append("\"minDiagonal\":").append(minDiagonal(symmetrizedSource)).append(",\n");
            sb.append("\"maxAbsDiagonal\":").append(maxAbsDiagonal(symmetrizedSource)).append(",\n");
            sb.append("\"gershgorinLowerBound\":").append(gershgorinLowerBound(symmetrizedSource)).append(",\n");
            sb.append("\"originalSource\":").append(jsonMatrix(originalSource)).append(",\n");
            sb.append("\"symmetrizedSource\":").append(jsonMatrix(symmetrizedSource));
            if (failureDump != null) {
                failureDump.appendJson(sb, context, originalSource, symmetrizedSource, jitterBase);
            } else {
                sb.append("\n");
            }
            sb.append("}\n");

            Files.write(Paths.get(path), sb.toString().getBytes(StandardCharsets.UTF_8));
            System.err.println("OU_SPD_FAILURE_DUMP " + path);
        } catch (final Exception e) {
            System.err.println("OU_SPD_FAILURE_DUMP_FAILED context=" + context + " reason=" + e.getMessage());
        }
    }

    private static void copySquare(final double[][] source, final double[][] destination, final int dimension) {
        for (int i = 0; i < dimension; ++i) {
            System.arraycopy(source[i], 0, destination[i], 0, dimension);
        }
    }

    private static void addDiagonal(final double[][] matrix, final int dimension, final double value) {
        for (int i = 0; i < dimension; ++i) {
            matrix[i][i] += value;
        }
    }

    public static void symmetrizeInPlace(final DenseMatrix64F matrix) {
        final int d = matrix.numRows;
        for (int i = 0; i < d; ++i) {
            for (int j = i + 1; j < d; ++j) {
                final double value = 0.5 * (matrix.unsafe_get(i, j) + matrix.unsafe_get(j, i));
                matrix.unsafe_set(i, j, value);
                matrix.unsafe_set(j, i, value);
            }
        }
    }

    public static void symmetrizeCopy(final DenseMatrix64F source,
                                      final DenseMatrix64F destination) {
        final int d = source.numRows;
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                destination.unsafe_set(i, j,
                        0.5 * (source.unsafe_get(i, j) + source.unsafe_get(j, i)));
            }
        }
    }

    public static void addDiagonalJitter(final DenseMatrix64F matrix,
                                         final double jitter) {
        for (int i = 0; i < matrix.numRows; ++i) {
            matrix.unsafe_set(i, i, matrix.unsafe_get(i, i) + jitter);
        }
    }

    private static void symmetrizeInPlace(final double[][] matrix, final int dimension) {
        for (int i = 0; i < dimension; ++i) {
            for (int j = i + 1; j < dimension; ++j) {
                final double value = 0.5 * (matrix[i][j] + matrix[j][i]);
                matrix[i][j] = value;
                matrix[j][i] = value;
            }
        }
    }

    private static double maxAbsDiagonal(final double[][] matrix, final int dimension) {
        double max = 0.0;
        for (int i = 0; i < dimension; ++i) {
            max = Math.max(max, Math.abs(matrix[i][i]));
        }
        return max;
    }

    public static double maxAbsDiagonal(final DenseMatrix64F matrix) {
        double max = 0.0;
        for (int i = 0; i < matrix.numRows; ++i) {
            max = Math.max(max, Math.abs(matrix.unsafe_get(i, i)));
        }
        return max;
    }

    private static double minDiagonal(final double[][] matrix, final int dimension) {
        double min = Double.POSITIVE_INFINITY;
        for (int i = 0; i < dimension; ++i) {
            min = Math.min(min, matrix[i][i]);
        }
        return min;
    }

    public static double minDiagonal(final DenseMatrix64F matrix) {
        final int d = Math.min(matrix.numRows, matrix.numCols);
        double min = Double.POSITIVE_INFINITY;
        for (int i = 0; i < d; ++i) {
            min = Math.min(min, matrix.unsafe_get(i, i));
        }
        return min;
    }

    private static double gershgorinLowerBound(final double[][] matrix, final int dimension) {
        double lowerBound = Double.POSITIVE_INFINITY;
        for (int i = 0; i < dimension; ++i) {
            double radius = 0.0;
            for (int j = 0; j < dimension; ++j) {
                if (i != j) {
                    radius += Math.abs(matrix[i][j]);
                }
            }
            lowerBound = Math.min(lowerBound, matrix[i][i] - radius);
        }
        return lowerBound;
    }

    public static double gershgorinLowerBound(final DenseMatrix64F matrix) {
        final int d = Math.min(matrix.numRows, matrix.numCols);
        double lowerBound = Double.POSITIVE_INFINITY;
        for (int i = 0; i < d; ++i) {
            double radius = 0.0;
            for (int j = 0; j < d; ++j) {
                if (i != j) {
                    radius += Math.abs(matrix.unsafe_get(i, j));
                }
            }
            lowerBound = Math.min(lowerBound, matrix.unsafe_get(i, i) - radius);
        }
        return lowerBound;
    }

    public static boolean isFinite(final double[] values) {
        for (double value : values) {
            if (!Double.isFinite(value)) {
                return false;
            }
        }
        return true;
    }

    public static boolean isFinite(final double[][] values) {
        for (double[] row : values) {
            if (!isFinite(row)) {
                return false;
            }
        }
        return true;
    }

    public static boolean isFinite(final DenseMatrix64F matrix) {
        final double[] data = matrix.getData();
        for (double value : data) {
            if (!Double.isFinite(value)) {
                return false;
            }
        }
        return true;
    }

    public static boolean hasNaN(final DenseMatrix64F matrix) {
        final double[] data = matrix.getData();
        for (double value : data) {
            if (Double.isNaN(value)) {
                return true;
            }
        }
        return false;
    }

    public static boolean hasInfinity(final DenseMatrix64F matrix) {
        final double[] data = matrix.getData();
        for (double value : data) {
            if (Double.isInfinite(value)) {
                return true;
            }
        }
        return false;
    }

    public static String summarizeDenseMatrix(final DenseMatrix64F matrix) {
        int nanCount = 0;
        int posInfCount = 0;
        int negInfCount = 0;
        double minFinite = Double.POSITIVE_INFINITY;
        double maxFinite = Double.NEGATIVE_INFINITY;
        final double[] data = matrix.getData();
        for (double value : data) {
            if (Double.isNaN(value)) {
                nanCount++;
            } else if (value == Double.POSITIVE_INFINITY) {
                posInfCount++;
            } else if (value == Double.NEGATIVE_INFINITY) {
                negInfCount++;
            } else {
                minFinite = Math.min(minFinite, value);
                maxFinite = Math.max(maxFinite, value);
            }
        }
        if (minFinite == Double.POSITIVE_INFINITY) {
            minFinite = Double.NaN;
            maxFinite = Double.NaN;
        }
        return "{rows=" + matrix.numRows
                + ",cols=" + matrix.numCols
                + ",nan=" + nanCount
                + ",posInf=" + posInfCount
                + ",negInf=" + negInfCount
                + ",minFinite=" + minFinite
                + ",maxFinite=" + maxFinite
                + "}";
    }

    public static String jsonMatrix(final DenseMatrix64F matrix) {
        final StringBuilder sb = new StringBuilder();
        sb.append("[");
        for (int i = 0; i < matrix.numRows; ++i) {
            if (i > 0) {
                sb.append(",");
            }
            sb.append("[");
            for (int j = 0; j < matrix.numCols; ++j) {
                if (j > 0) {
                    sb.append(",");
                }
                sb.append(matrix.unsafe_get(i, j));
            }
            sb.append("]");
        }
        sb.append("]");
        return sb.toString();
    }

    public static String jsonVector(final DenseMatrix64F vector) {
        final StringBuilder sb = new StringBuilder();
        sb.append("[");
        for (int i = 0; i < vector.numRows; ++i) {
            if (i > 0) {
                sb.append(",");
            }
            sb.append(vector.unsafe_get(i, 0));
        }
        sb.append("]");
        return sb.toString();
    }

    private static boolean isFinite(final double[][] matrix, final int dimension) {
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                if (!Double.isFinite(matrix[i][j])) {
                    return false;
                }
            }
        }
        return true;
    }

    private static void checkSquare(final double[][] matrix, final int dimension, final String label) {
        if (matrix.length < dimension) {
            throw new IllegalArgumentException(label + " has too few rows.");
        }
        for (int i = 0; i < dimension; ++i) {
            if (matrix[i].length < dimension) {
                throw new IllegalArgumentException(label + "[" + i + "] has too few columns.");
            }
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
