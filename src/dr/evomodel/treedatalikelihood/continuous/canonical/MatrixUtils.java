/*
 * MatrixUtils.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 *
 */

package dr.evomodel.treedatalikelihood.continuous.canonical;

import dr.evomodel.treedatalikelihood.continuous.canonical.math.MatrixOps;
import dr.math.matrixAlgebra.missingData.InversionResult;
import dr.math.matrixAlgebra.missingData.MissingOps;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import java.util.Arrays;

/**
 * Stateless matrix utilities used internally by the framework.
 *
 * <p>All methods operate on row-major {@code double[]} arrays.
 *
 * <p>Most methods in this class are superseded by
 * {@link dr.evomodel.treedatalikelihood.continuous.canonical.math.MatrixOps}, which
 * provides a cleaner, unified API. New code should use {@code MatrixOps} directly.
 * {@link #invertSymmetricPositiveDefiniteCompact} is retained here because it has
 * unique semantics (in-place Cholesky factorization of the source buffer, explicit
 * work arrays) that are not replicated in {@code MatrixOps}.
 *
 * @author Filippo Monti
 * @author Marc A. Suchard
 */
public final class MatrixUtils {

    private MatrixUtils() { }

    // -----------------------------------------------------------------------
    // Inversion
    // -----------------------------------------------------------------------

    /**
     * Inverts the symmetric positive-definite matrix {@code src} and writes the result into
     * {@code dst}. Both arrays are row-major, length {@code dim × dim}.
     *
     * @deprecated Use {@link MatrixOps#invertSPD(double[], double[], int)} instead.
     */
    @Deprecated
    public static void invertSymmetric(double[] src, double[] dst, int dim) {
        MatrixOps.invertSPD(src, dst, dim);
    }

    /**
     * Inverts a compact symmetric positive-definite matrix without dynamic allocation.
     *
     * <p>The leading {@code dim x dim} block of {@code src} is interpreted as a row-major
     * compact matrix and is overwritten with its lower-triangular Cholesky factor. The
     * leading {@code dim x dim} block of {@code dst} receives the inverse matrix.
     *
     * <p>The returned value is {@code log |src|}, i.e. the log-determinant of the original
     * variance matrix before overwrite.
     *
     * @param src row-major compact SPD matrix, overwritten in-place with its Cholesky factor
     * @param dst output inverse matrix buffer; may be larger than {@code dim x dim}
     * @param dim compact matrix dimension
     * @param forwardWork temporary vector of length at least {@code dim}
     * @param backwardWork temporary vector of length at least {@code dim}
     * @return {@code log |src|}
     */
    public static double invertSymmetricPositiveDefiniteCompact(final double[] src,
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

        for (int row = 0; row < dim; ++row) {
            final int rowOffset = row * dim;
            for (int col = row + 1; col < dim; ++col) {
                final double sym = 0.5 * (dst[rowOffset + col] + dst[col * dim + row]);
                dst[rowOffset + col] = sym;
                dst[col * dim + row] = sym;
            }
        }

        return logDeterminant;
    }

    /**
     * Inverts a symmetric variance matrix that may contain {@code +∞} entries on the diagonal
     * (representing missing trait dimensions), and returns an {@link InversionResult} with the
     * log-determinant of the resulting <em>precision</em> matrix and the effective dimension.
     *
     * <p>Missing dimensions (∞ diagonal in variance) are mapped to 0 precision; observed
     * dimensions are inverted in the usual way. Off-diagonal entries for missing rows/columns
     * are set to zero in the output.
     *
     * @param src variance matrix (row-major, length {@code dim×dim}); may have ∞ entries
     * @param dst output precision matrix (row-major, length {@code dim×dim})
     * @param dim matrix dimension
     * @return inversion result with effective dimension and log-determinant of the precision
     */
    /** @deprecated Use {@link MatrixOps#safeInvertVariance(double[], double[], int)} instead. */
    @Deprecated
    public static InversionResult safeInvertVariance(double[] src, double[] dst, int dim) {
        return MatrixOps.safeInvertVariance(src, dst, dim);
    }

    /**
     * Inverts a symmetric precision matrix that may contain {@code 0} entries on the diagonal
     * (representing missing trait dimensions), and returns an {@link InversionResult} with the
     * log-determinant of the resulting <em>precision</em> (not variance) and the effective dimension.
     *
     * <p>Missing dimensions (0 diagonal in precision) are mapped to ∞ variance; observed
     * dimensions are inverted in the usual way.
     *
     * @param src precision matrix (row-major, length {@code dim×dim}); may have 0 entries
     * @param dst output variance matrix (row-major, length {@code dim×dim})
     * @param dim matrix dimension
     * @return inversion result with effective dimension and log-determinant of the <em>precision</em>
     */
    /** @deprecated Use {@link MatrixOps#safeInvertPrecision(double[], double[], int)} instead. */
    @Deprecated
    public static InversionResult safeInvertPrecision(double[] src, double[] dst, int dim) {
        return MatrixOps.safeInvertPrecision(src, dst, dim);
    }

    // -----------------------------------------------------------------------
    // Log-determinant
    // -----------------------------------------------------------------------

    /**
     * Returns the log-determinant of a symmetric positive-definite matrix via Cholesky
     * decomposition.
     *
     * @param matrix row-major symmetric PD matrix, length {@code dim × dim}
     * @param dim    matrix dimension
     * @return {@code log |matrix|}
     */
    /** @deprecated Use {@link MatrixOps#logDeterminant(double[], int)} instead. */
    @Deprecated
    public static double logDeterminant(double[] matrix, int dim) {
        return MatrixOps.logDeterminant(matrix, dim);
    }

    // -----------------------------------------------------------------------
    // Matrix–vector product
    // -----------------------------------------------------------------------

    /**
     * Computes the matrix–vector product {@code y = A · x} and stores the result in {@code y}.
     *
     * @param A   row-major matrix, length {@code m × n}
     * @param x   input vector, length {@code n}
     * @param y   output vector, length {@code m} (must be pre-allocated)
     * @param m   number of rows
     * @param n   number of columns
     */
    /** @deprecated Use {@link MatrixOps#matVec(double[], double[], double[], int, int)} instead. */
    @Deprecated
    public static void matVec(double[] A, double[] x, double[] y, int m, int n) {
        MatrixOps.matVec(A, x, y, m, n);
    }

    // -----------------------------------------------------------------------
    // Matrix–matrix product
    // -----------------------------------------------------------------------

    /**
     * Computes {@code C = A · B} (row-major) and stores the result in {@code C}.
     *
     * @param A   row-major matrix, length {@code m × k}
     * @param B   row-major matrix, length {@code k × n}
     * @param C   output matrix, length {@code m × n} (must be pre-allocated, may not alias A or B)
     * @param m   rows of A (and C)
     * @param k   cols of A / rows of B
     * @param n   cols of B (and C)
     */
    /** @deprecated Use {@link MatrixOps#matMul(double[], double[], double[], int, int, int)} instead. */
    @Deprecated
    public static void matMul(double[] A, double[] B, double[] C, int m, int k, int n) {
        MatrixOps.matMul(A, B, C, m, k, n);
    }

    // -----------------------------------------------------------------------
    // Scale / add
    // -----------------------------------------------------------------------

    /**
     * Computes {@code dst = src + scalar * add} element-wise, storing the result in {@code dst}.
     * All arrays must have the same length.
     */
    /** @deprecated Use {@link MatrixOps#addScaled(double[], double[], double[], double, int)} instead. */
    @Deprecated
    public static void addScaled(double[] dst, double[] src, double[] add, double scalar, int len) {
        MatrixOps.addScaled(dst, src, add, scalar, len);
    }

    /** @deprecated Use {@link MatrixOps#scaleInPlace(double[], double, int)} instead. */
    @Deprecated
    public static void scaleInPlace(double[] arr, double scalar, int len) {
        MatrixOps.scaleInPlace(arr, scalar, len);
    }

    // -----------------------------------------------------------------------
    // Quadratic form
    // -----------------------------------------------------------------------

    /**
     * Computes the quadratic form {@code x^T · A · x}.
     *
     * @param x   vector, length {@code d}
     * @param A   row-major matrix, length {@code d × d}
     * @param d   dimension
     * @param tmp scratch vector, length {@code d}
     * @return {@code x^T A x}
     */
    /** @deprecated Use {@link MatrixOps#quadraticForm(double[], double[], int, double[])} instead. */
    @Deprecated
    public static double quadraticForm(double[] x, double[] A, int d, double[] tmp) {
        return MatrixOps.quadraticForm(x, A, d, tmp);
    }

    // -----------------------------------------------------------------------
    // Difference
    // -----------------------------------------------------------------------

    /**
     * Fills {@code diff} with {@code a - b} element-wise.
     */
    /** @deprecated Use {@link MatrixOps#subtract(double[], double[], double[], int)} instead. */
    @Deprecated
    public static void subtract(double[] diff, double[] a, double[] b, int len) {
        MatrixOps.subtract(diff, a, b, len);
    }
}
