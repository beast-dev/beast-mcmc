package dr.evomodel.treedatalikelihood.continuous.gaussian.message;

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

    private static double minDiagonal(final double[][] matrix, final int dimension) {
        double min = Double.POSITIVE_INFINITY;
        for (int i = 0; i < dimension; ++i) {
            min = Math.min(min, matrix[i][i]);
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
}
