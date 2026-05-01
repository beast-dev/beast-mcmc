package dr.evomodel.treedatalikelihood.continuous.canonical.math;

import java.util.Arrays;

/**
 * Stateless flat row-major matrix operations for the canonical diffusion framework.
 *
 * <p>All methods operate on {@code double[]} arrays in row-major order:
 * element {@code (i,j)} of an {@code m×n} matrix is at index {@code i*n + j}.
 * No {@code double[][]} appears in any method signature.
 *
 * <p>This class supersedes
 * {@link dr.evomodel.treedatalikelihood.continuous.canonical.MatrixUtils} and the
 * flat-array portion of
 * {@link dr.evomodel.treedatalikelihood.continuous.canonical.message.GaussianMatrixOps}.
 * Those classes retain their implementations for backward compatibility but should not
 * be used in new code.
 *
 * <p>Conversion helpers ({@code toFlat}, {@code fromFlat}, …) are provided solely for
 * use at {@code double[][]} boundaries (external BEAST model parameters). Do not store
 * matrices as {@code double[][]} inside the canonical framework.
 */
public final class MatrixOps {

    public static final double LOG_TWO_PI = Math.log(2.0 * Math.PI);
    public static final double MIN_DIAGONAL_JITTER = 1e-10;

    private MatrixOps() { }

    // -----------------------------------------------------------------------
    // Matrix–vector  (y = A x)
    // -----------------------------------------------------------------------

    /** {@code y = A x}, A is {@code m×n} row-major. */
    public static void matVec(double[] A, double[] x, double[] y, int m, int n) {
        for (int i = 0; i < m; i++) {
            double s = 0.0;
            final int base = i * n;
            for (int j = 0; j < n; j++) s += A[base + j] * x[j];
            y[i] = s;
        }
    }

    /** {@code y = A x}, A is {@code dim×dim} square. */
    public static void matVec(double[] A, double[] x, double[] y, int dim) {
        matVec(A, x, y, dim, dim);
    }

    // -----------------------------------------------------------------------
    // Matrix–matrix  (C = A B)
    // -----------------------------------------------------------------------

    /**
     * {@code C = A B}, A is {@code m×k}, B is {@code k×n}, C is {@code m×n}, all row-major.
     * C must not alias A or B.
     */
    public static void matMul(double[] A, double[] B, double[] C, int m, int k, int n) {
        for (int i = 0; i < m; i++) {
            final int rA = i * k, rC = i * n;
            for (int j = 0; j < n; j++) {
                double s = 0.0;
                for (int l = 0; l < k; l++) s += A[rA + l] * B[l * n + j];
                C[rC + j] = s;
            }
        }
    }

    /** {@code C = A B}, all {@code dim×dim} square. */
    public static void matMul(double[] A, double[] B, double[] C, int dim) {
        matMul(A, B, C, dim, dim, dim);
    }

    /** {@code C = A B^T}, all {@code dim×dim} square. */
    public static void matMulTransposedRight(double[] A, double[] B, double[] C, int dim) {
        for (int i = 0; i < dim; i++) {
            final int iOff = i * dim;
            for (int j = 0; j < dim; j++) {
                double sum = 0.0;
                final int jOff = j * dim;
                for (int k = 0; k < dim; k++) {
                    sum += A[iOff + k] * B[jOff + k];
                }
                C[iOff + j] = sum;
            }
        }
    }

    // -----------------------------------------------------------------------
    // Misc flat arithmetic
    // -----------------------------------------------------------------------

    /** Symmetrize in place: {@code m[i,j] = m[j,i] = 0.5*(m[i,j]+m[j,i])}. */
    public static void symmetrize(double[] m, int dim) {
        for (int i = 0; i < dim; i++) {
            for (int j = i + 1; j < dim; j++) {
                final double avg = 0.5 * (m[i * dim + j] + m[j * dim + i]);
                m[i * dim + j] = avg;
                m[j * dim + i] = avg;
            }
        }
    }

    /** {@code acc[i] += inc[i]} for {@code i} in {@code [0, len)}. */
    public static void addInPlace(double[] acc, double[] inc, int len) {
        for (int i = 0; i < len; i++) acc[i] += inc[i];
    }

    /** {@code dst[i] = src[i] + scalar * add[i]}. */
    public static void addScaled(double[] dst, double[] src, double[] add, double scalar, int len) {
        for (int i = 0; i < len; i++) dst[i] = src[i] + scalar * add[i];
    }

    /** {@code arr[i] *= scalar}. */
    public static void scaleInPlace(double[] arr, double scalar, int len) {
        for (int i = 0; i < len; i++) arr[i] *= scalar;
    }

    /** {@code diff[i] = a[i] - b[i]}. */
    public static void subtract(double[] diff, double[] a, double[] b, int len) {
        for (int i = 0; i < len; i++) diff[i] = a[i] - b[i];
    }

    /** {@code out[i,j] = left[i] * right[j]}. */
    public static void outerProduct(double[] left, double[] right, double[] out, int dim) {
        for (int i = 0; i < dim; i++) {
            final int row = i * dim;
            for (int j = 0; j < dim; j++) {
                out[row + j] = left[i] * right[j];
            }
        }
    }

    /**
     * Quadratic form {@code x^T A x}.
     *
     * @param tmp scratch vector of length {@code dim}; overwritten on return
     */
    public static double quadraticForm(double[] x, double[] A, int dim, double[] tmp) {
        matVec(A, x, tmp, dim);
        double r = 0.0;
        for (int i = 0; i < dim; i++) r += x[i] * tmp[i];
        return r;
    }

    /** {@code System.arraycopy(src, 0, dst, 0, dim*dim)}. */
    public static void copyMatrix(double[] src, double[] dst, int dim) {
        System.arraycopy(src, 0, dst, 0, dim * dim);
    }

    /** {@code System.arraycopy(src, 0, dst, 0, src.length)}. */
    public static void copyVector(double[] src, double[] dst) {
        System.arraycopy(src, 0, dst, 0, src.length);
    }

    /** Transpose a square row-major matrix in place. */
    public static void transposeInPlace(double[] m, int dim) {
        for (int i = 0; i < dim; i++) {
            for (int j = i + 1; j < dim; j++) {
                final int ij = i * dim + j, ji = j * dim + i;
                final double tmp = m[ij];
                m[ij] = m[ji];
                m[ji] = tmp;
            }
        }
    }

    /**
     * Transpose: {@code dst[j*dim+i] = src[i*dim+j]}.
     * {@code dst} must not alias {@code src}.
     */
    public static void transpose(double[] src, double[] dst, int dim) {
        for (int i = 0; i < dim; i++)
            for (int j = 0; j < dim; j++)
                dst[j * dim + i] = src[i * dim + j];
    }

    // -----------------------------------------------------------------------
    // Cholesky factorization
    // -----------------------------------------------------------------------

    /**
     * Attempts to compute the lower Cholesky factor of a symmetric positive-definite
     * row-major matrix {@code matrix} into {@code lowerOut}.
     * {@code matrix} is not modified.
     *
     * @return {@code true} on success; {@code false} if the matrix is not positive definite
     */
    public static boolean tryCholesky(double[] matrix, double[] lowerOut, int dim) {
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j <= i; j++) {
                double s = matrix[i * dim + j];
                for (int k = 0; k < j; k++) s -= lowerOut[i * dim + k] * lowerOut[j * dim + k];
                if (i == j) {
                    if (s <= 0.0) return false;
                    lowerOut[i * dim + j] = Math.sqrt(s);
                } else {
                    final double d = lowerOut[j * dim + j];
                    if (d == 0.0) return false;
                    lowerOut[i * dim + j] = s / d;
                }
            }
            for (int j = i + 1; j < dim; j++) lowerOut[i * dim + j] = 0.0;
        }
        return true;
    }

    /**
     * Given a lower Cholesky factor {@code lower} (row-major, {@code dim×dim}),
     * computes {@code L^{-1}} into {@code lowerInvOut} and the SPD inverse
     * {@code A^{-1} = L^{-T} L^{-1}} into {@code out}.
     *
     * @param lower      Cholesky lower factor L, row-major
     * @param lowerInvOut scratch for L^{-1}, row-major; must be length {@code dim*dim}
     * @param out        output A^{-1}, row-major; must be length {@code dim*dim}
     * @param dim        matrix dimension
     * @return log-determinant of A (= 2 * sum of log-diagonals of L)
     */
    public static double invertFromCholesky(double[] lower, double[] lowerInvOut, double[] out, int dim) {
        double logDet = 0.0;
        for (int i = 0; i < dim; i++) logDet += Math.log(lower[i * dim + i]);
        logDet *= 2.0;

        for (int col = 0; col < dim; col++) {
            for (int row = 0; row < dim; row++) {
                double s = (row == col) ? 1.0 : 0.0;
                for (int k = 0; k < row; k++) s -= lower[row * dim + k] * lowerInvOut[k * dim + col];
                lowerInvOut[row * dim + col] = s / lower[row * dim + row];
            }
        }

        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                double s = 0.0;
                for (int k = 0; k < dim; k++) s += lowerInvOut[k * dim + i] * lowerInvOut[k * dim + j];
                out[i * dim + j] = s;
            }
        }
        symmetrize(out, dim);
        return logDet;
    }

    /**
     * Inverts the leading compact {@code dim x dim} row-major SPD block of {@code src}
     * into {@code dst}. The leading block of {@code src} is overwritten with the
     * Cholesky factor. This allocation-free variant exists for reduced missing-data
     * matrices whose active dimension changes inside a larger scratch buffer.
     *
     * @return log determinant of the original compact SPD block
     */
    public static double invertSPDCompact(final double[] src,
                                          final double[] dst,
                                          final int dim,
                                          final double[] forwardWork,
                                          final double[] backwardWork) {
        if (forwardWork.length < dim || backwardWork.length < dim) {
            throw new IllegalArgumentException("workspace vectors are too small for compact inversion");
        }
        if (src.length < dim * dim || dst.length < dim * dim) {
            throw new IllegalArgumentException("matrix buffers are too small for compact inversion");
        }

        double logDeterminant = 0.0;
        for (int row = 0; row < dim; ++row) {
            final int rowOffset = row * dim;
            for (int col = 0; col <= row; ++col) {
                double sum = src[rowOffset + col];
                final int colOffset = col * dim;
                for (int k = 0; k < col; ++k) {
                    sum -= src[rowOffset + k] * src[colOffset + k];
                }
                if (row == col) {
                    if (!(sum > 0.0) || Double.isNaN(sum)) {
                        throw new IllegalArgumentException("Matrix is not symmetric positive definite.");
                    }
                    final double diagonal = Math.sqrt(sum);
                    src[rowOffset + row] = diagonal;
                    logDeterminant += 2.0 * Math.log(diagonal);
                } else {
                    src[rowOffset + col] = sum / src[colOffset + col];
                }
            }
            Arrays.fill(src, rowOffset + row + 1, rowOffset + dim, 0.0);
        }

        Arrays.fill(dst, 0, dim * dim, 0.0);
        for (int column = 0; column < dim; ++column) {
            for (int row = 0; row < dim; ++row) {
                double sum = row == column ? 1.0 : 0.0;
                final int rowOffset = row * dim;
                for (int k = 0; k < row; ++k) {
                    sum -= src[rowOffset + k] * forwardWork[k];
                }
                forwardWork[row] = sum / src[rowOffset + row];
            }

            for (int row = dim - 1; row >= 0; --row) {
                double sum = forwardWork[row];
                for (int k = row + 1; k < dim; ++k) {
                    sum -= src[k * dim + row] * backwardWork[k];
                }
                backwardWork[row] = sum / src[row * dim + row];
            }

            for (int row = 0; row < dim; ++row) {
                dst[row * dim + column] = backwardWork[row];
            }
        }

        symmetrize(dst, dim);
        return logDeterminant;
    }

    // -----------------------------------------------------------------------
    // SPD inversion and missing-aware safe variants
    // -----------------------------------------------------------------------

    /**
     * Invert a symmetric positive-definite matrix using flat Cholesky workspaces.
     * Prefer {@link GaussianFormConverter} for inversion in the canonical hot path.
     *
     * @throws IllegalArgumentException if the matrix is not symmetric positive definite
     */
    public static void invertSPD(double[] src, double[] dst, int dim) {
        final double[] lower = new double[dim * dim];
        final double[] lowerInverse = new double[dim * dim];
        if (!tryCholesky(src, lower, dim)) {
            throw new IllegalArgumentException("Matrix is not symmetric positive definite.");
        }
        invertFromCholesky(lower, lowerInverse, dst, dim);
    }

    /**
     * Invert a variance matrix that may contain ∞ on the diagonal (missing traits),
     * returning the precision matrix and a {@link CanonicalInversionResult} with the
     * effective dimension and log-determinant of the finite variance block.
     */
    public static CanonicalInversionResult safeInvertVariance(double[] src, double[] dst, int dim) {
        return safeInvertWithMissing(src, dst, dim, true);
    }

    /**
     * Invert a precision matrix that may contain 0 on the diagonal (missing traits),
     * returning the variance matrix and a {@link CanonicalInversionResult} with the
     * effective dimension and log-determinant of the finite precision block.
     */
    public static CanonicalInversionResult safeInvertPrecision(double[] src, double[] dst, int dim) {
        return safeInvertWithMissing(src, dst, dim, false);
    }

    /**
     * Log-determinant via flat Cholesky decomposition.
     */
    public static double logDeterminant(double[] m, int dim) {
        final double[] lower = new double[dim * dim];
        if (!tryCholesky(m, lower, dim)) {
            throw new IllegalArgumentException("Matrix is not symmetric positive definite.");
        }
        double logDet = 0.0;
        for (int i = 0; i < dim; i++) {
            logDet += Math.log(lower[i * dim + i]);
        }
        return 2.0 * logDet;
    }

    private static CanonicalInversionResult safeInvertWithMissing(final double[] src,
                                                                  final double[] dst,
                                                                  final int dim,
                                                                  final boolean inputIsVariance) {
        Arrays.fill(dst, 0, dim * dim, 0.0);

        final int[] finiteIndices = new int[dim];
        int finiteCount = 0;
        int exactCount = 0;

        for (int i = 0; i < dim; i++) {
            final double diagonal = src[i * dim + i];
            if (Double.isNaN(diagonal) || diagonal == Double.NEGATIVE_INFINITY) {
                throw new IllegalArgumentException("Matrix diagonal contains an unsupported missing-data marker.");
            }

            if (isFiniteNonZero(diagonal)) {
                finiteIndices[finiteCount++] = i;
            } else if (diagonal == 0.0) {
                dst[i * dim + i] = Double.POSITIVE_INFINITY;
                if (inputIsVariance) {
                    exactCount++;
                }
            } else if (diagonal == Double.POSITIVE_INFINITY) {
                if (!inputIsVariance) {
                    exactCount++;
                }
            } else {
                throw new IllegalArgumentException("Matrix diagonal contains an unsupported missing-data marker.");
            }
        }

        final int effectiveDimension = finiteCount + exactCount;
        double logDeterminant;

        if (finiteCount > 0) {
            final double[] compactSource = new double[finiteCount * finiteCount];
            final double[] compactInverse = new double[finiteCount * finiteCount];
            final double[] forwardWork = new double[finiteCount];
            final double[] backwardWork = new double[finiteCount];

            for (int row = 0; row < finiteCount; row++) {
                final int sourceRow = finiteIndices[row];
                for (int col = 0; col < finiteCount; col++) {
                    compactSource[row * finiteCount + col] = src[sourceRow * dim + finiteIndices[col]];
                }
            }

            final double sourceLogDeterminant = invertSPDCompact(
                    compactSource, compactInverse, finiteCount, forwardWork, backwardWork);
            for (int row = 0; row < finiteCount; row++) {
                final int targetRow = finiteIndices[row];
                for (int col = 0; col < finiteCount; col++) {
                    dst[targetRow * dim + finiteIndices[col]] = compactInverse[row * finiteCount + col];
                }
            }

            logDeterminant = sourceLogDeterminant;
        } else {
            logDeterminant = Double.NEGATIVE_INFINITY;
        }

        if (exactCount > 0) {
            logDeterminant = Double.POSITIVE_INFINITY;
        }

        return new CanonicalInversionResult(
                CanonicalInversionResult.getCode(dim, effectiveDimension),
                effectiveDimension,
                logDeterminant);
    }

    private static boolean isFiniteNonZero(final double value) {
        return value != 0.0 && !Double.isNaN(value) && !Double.isInfinite(value);
    }

    // -----------------------------------------------------------------------
    // Conversion helpers — double[][] boundary only
    //
    // Use these only when receiving data from external BEAST model objects that
    // provide double[][]. Do not introduce new double[][] storage in canonical code.
    // -----------------------------------------------------------------------

    /** Copy {@code double[][]} into row-major {@code double[]}. */
    public static void toFlat(double[][] src, double[] dst, int dim) {
        for (int i = 0; i < dim; i++) System.arraycopy(src[i], 0, dst, i * dim, dim);
    }

    /** Copy row-major {@code double[]} into {@code double[][]}. */
    public static void fromFlat(double[] src, double[][] dst, int dim) {
        for (int i = 0; i < dim; i++) System.arraycopy(src, i * dim, dst[i], 0, dim);
    }

    /**
     * Copy the transpose of a row-major {@code double[]} into a {@code double[][]}.
     * i.e. {@code dst[j][i] = src[i*dim+j]}.
     */
    public static void transposedFromFlat(double[] src, double[][] dst, int dim) {
        for (int i = 0; i < dim; i++)
            for (int j = 0; j < dim; j++)
                dst[j][i] = src[i * dim + j];
    }

    /** Copy a {@code double[][]} row-major into a column-major parameter vector. */
    public static void toColumnMajorParam(double[][] src, double[] dst, int dim) {
        int idx = 0;
        for (int col = 0; col < dim; col++)
            for (int row = 0; row < dim; row++)
                dst[idx++] = src[row][col];
    }

    /** Copy a row-major flat matrix into a column-major parameter vector. */
    public static void toColumnMajorParam(double[] src, double[] dst, int dim) {
        for (int row = 0; row < dim; row++)
            for (int col = 0; col < dim; col++)
                dst[col * dim + row] = src[row * dim + col];
    }

    /** Copy a column-major parameter vector into a row-major flat matrix. */
    public static void fromColumnMajorParam(double[] src, double[] dst, int dim) {
        for (int row = 0; row < dim; row++)
            for (int col = 0; col < dim; col++)
                dst[row * dim + col] = src[col * dim + row];
    }
}
