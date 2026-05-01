package dr.evomodel.treedatalikelihood.continuous.canonical.math;

import dr.math.matrixAlgebra.missingData.InversionResult;
import dr.math.matrixAlgebra.missingData.MissingOps;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

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

    // -----------------------------------------------------------------------
    // EJML-backed inversion (for BEAST parameter boundaries and safe variants)
    // -----------------------------------------------------------------------

    /**
     * Invert a symmetric positive-definite matrix using EJML's Cholesky solver.
     * Prefer {@link GaussianFormConverter} for inversion in the canonical hot path.
     *
     * @throws IllegalArgumentException if the matrix is singular or ill-conditioned
     */
    public static void invertSPD(double[] src, double[] dst, int dim) {
        final DenseMatrix64F m = DenseMatrix64F.wrap(dim, dim, src);
        final DenseMatrix64F inv = new DenseMatrix64F(dim, dim);
        if (!CommonOps.invert(m, inv))
            throw new IllegalArgumentException("Matrix is not invertible (singular or ill-conditioned).");
        System.arraycopy(inv.data, 0, dst, 0, dim * dim);
    }

    /**
     * Invert a variance matrix that may contain ∞ on the diagonal (missing traits),
     * returning the precision matrix and an {@link InversionResult} with the
     * effective dimension and log-determinant of the precision.
     */
    public static InversionResult safeInvertVariance(double[] src, double[] dst, int dim) {
        final DenseMatrix64F srcM = new DenseMatrix64F(dim, dim);
        System.arraycopy(src, 0, srcM.data, 0, dim * dim);
        final DenseMatrix64F dstM = new DenseMatrix64F(dim, dim);
        final InversionResult result = MissingOps.safeInvertVariance(srcM, dstM, true);
        System.arraycopy(dstM.data, 0, dst, 0, dim * dim);
        return result;
    }

    /**
     * Invert a precision matrix that may contain 0 on the diagonal (missing traits),
     * returning the variance matrix and an {@link InversionResult} with the
     * effective dimension and log-determinant of the original precision.
     */
    public static InversionResult safeInvertPrecision(double[] src, double[] dst, int dim) {
        final DenseMatrix64F srcM = new DenseMatrix64F(dim, dim);
        System.arraycopy(src, 0, srcM.data, 0, dim * dim);
        final DenseMatrix64F dstM = new DenseMatrix64F(dim, dim);
        final InversionResult result = MissingOps.safeInvertPrecision(srcM, dstM, true);
        System.arraycopy(dstM.data, 0, dst, 0, dim * dim);
        return result;
    }

    /**
     * Log-determinant via LU decomposition.
     * Use Cholesky-based determinant via {@link #invertFromCholesky} for SPD matrices.
     */
    public static double logDeterminant(double[] m, int dim) {
        return Math.log(Math.abs(CommonOps.det(DenseMatrix64F.wrap(dim, dim, m))));
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
