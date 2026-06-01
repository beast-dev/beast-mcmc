package dr.evomodel.continuous.ou;

import java.util.Arrays;

/**
 * Dense row-major matrix helpers for exact OU transitions.
 *
 * <p>This utility is intentionally self-contained so exact OU transition logic can live
 * without adding a new linear-algebra dependency. The matrix exponential uses
 * scaling-and-squaring with the order-13 Pade approximant.
 *
 * <p>For selection-matrix gradients we expose an <em>adjoint</em> (backpropagation) helper
 * for the matrix exponential. If
 * <pre>
 *     Y = expm(X)
 * </pre>
 * and the upstream sensitivity is G = dL/dY, then the Frobenius-adjoint sensitivity wrt X is
 * <pre>
 *     dL/dX = L_exp(X^T, G)^T,
 * </pre>
 * where L_exp is the Frechet derivative of the matrix exponential. The implementation
 * evaluates this adjoint via the standard block-exponential identity rather than by explicit
 * forward-direction loops over basis matrices.
 */
public final class MatrixExponentialUtils {

    private static final double THETA_13 = 5.371920351148152;
    private static final double[] PADE_13 = {
            64764752532480000.0,
            32382376266240000.0,
            7771770303897600.0,
            1187353796428800.0,
            129060195264000.0,
            10559470521600.0,
            670442572800.0,
            33522128640.0,
            1323241920.0,
            40840800.0,
            960960.0,
            16380.0,
            182.0,
            1.0
    };

    private MatrixExponentialUtils() {
        // no instances
    }

    public static void setIdentityFlat(final double[] matrix, final int dimension) {
        setIdentityFlat(matrix, 0, dimension);
    }

    public static void setIdentityFlat(final double[] matrix,
                                       final int offset,
                                       final int dimension) {
        Arrays.fill(matrix, offset, offset + dimension * dimension, 0.0);
        for (int i = 0; i < dimension; ++i) {
            matrix[offset + i * dimension + i] = 1.0;
        }
    }

    public static void fillFlat(final double[] matrix,
                                final int offset,
                                final int length,
                                final double value) {
        Arrays.fill(matrix, offset, offset + length, value);
    }

    public static void transposeFlat(final double[] source,
                                     final double[] out,
                                     final int dimension) {
        transposeFlat(source, 0, out, 0, dimension);
    }

    public static void transposeFlat(final double[] source,
                                     final int sourceOffset,
                                     final double[] out,
                                     final int outOffset,
                                     final int dimension) {
        for (int i = 0; i < dimension; ++i) {
            final int sourceRow = sourceOffset + i * dimension;
            for (int j = 0; j < dimension; ++j) {
                out[outOffset + j * dimension + i] = source[sourceRow + j];
            }
        }
    }

    public static void multiplyFlat(final double[] left,
                                    final double[] right,
                                    final double[] out,
                                    final int dimension) {
        multiplyFlat(left, 0, right, 0, out, 0, dimension);
    }

    public static void multiplyFlat(final double[] left,
                                    final int leftOffset,
                                    final double[] right,
                                    final int rightOffset,
                                    final double[] out,
                                    final int outOffset,
                                    final int dimension) {
        for (int i = 0; i < dimension; ++i) {
            final int leftRow = leftOffset + i * dimension;
            final int outRow = outOffset + i * dimension;
            for (int j = 0; j < dimension; ++j) {
                double sum = 0.0;
                for (int k = 0; k < dimension; ++k) {
                    sum += left[leftRow + k] * right[rightOffset + k * dimension + j];
                }
                out[outRow + j] = sum;
            }
        }
    }

    public static void symmetrizeFlat(final double[] matrix, final int dimension) {
        symmetrizeFlat(matrix, 0, dimension);
    }

    public static void symmetrizeFlat(final double[] matrix,
                                      final int offset,
                                      final int dimension) {
        for (int i = 0; i < dimension; ++i) {
            for (int j = i + 1; j < dimension; ++j) {
                final int ij = offset + i * dimension + j;
                final int ji = offset + j * dimension + i;
                final double average = 0.5 * (matrix[ij] + matrix[ji]);
                matrix[ij] = average;
                matrix[ji] = average;
            }
        }
    }

    public static double norm1Flat(final double[] matrix, final int dimension) {
        return norm1Flat(matrix, 0, dimension);
    }

    public static double norm1Flat(final double[] matrix,
                                   final int offset,
                                   final int dimension) {
        double best = 0.0;
        for (int j = 0; j < dimension; ++j) {
            double sum = 0.0;
            for (int i = 0; i < dimension; ++i) {
                sum += Math.abs(matrix[offset + i * dimension + j]);
            }
            if (sum > best) {
                best = sum;
            }
        }
        return best;
    }

    public static void expmFlat(final double[] matrix,
                                final double[] out,
                                final int dimension) {
        expmFlat(matrix, 0, out, 0, dimension);
    }

    public static void expmFlat(final double[] matrix,
                                final int matrixOffset,
                                final double[] out,
                                final int outOffset,
                                final int dimension) {
        final int length = dimension * dimension;
        final double[] a = new double[length];
        System.arraycopy(matrix, matrixOffset, a, 0, length);

        final double aNorm = norm1Flat(a, dimension);
        final int squarings = aNorm <= THETA_13 ? 0 : Math.max(0, (int) Math.ceil(log2(aNorm / THETA_13)));
        if (squarings > 0) {
            final double scale = Math.scalb(1.0, -squarings);
            for (int i = 0; i < length; ++i) {
                a[i] *= scale;
            }
        }

        final double[] a2 = new double[length];
        final double[] a4 = new double[length];
        final double[] a6 = new double[length];
        multiplyFlat(a, a, a2, dimension);
        multiplyFlat(a2, a2, a4, dimension);
        multiplyFlat(a2, a4, a6, dimension);

        final double[] tmp1 = new double[length];
        final double[] tmp2 = new double[length];
        final double[] uInner = new double[length];
        final double[] u = new double[length];
        final double[] v = new double[length];
        final double[] identity = new double[length];
        setIdentityFlat(identity, dimension);

        addScaledInPlaceFlat(tmp1, a6, PADE_13[13]);
        addScaledInPlaceFlat(tmp1, a4, PADE_13[11]);
        addScaledInPlaceFlat(tmp1, a2, PADE_13[9]);
        multiplyFlat(a6, tmp1, tmp2, dimension);
        addScaledInPlaceFlat(tmp2, a6, PADE_13[7]);
        addScaledInPlaceFlat(tmp2, a4, PADE_13[5]);
        addScaledInPlaceFlat(tmp2, a2, PADE_13[3]);
        addScaledInPlaceFlat(tmp2, identity, PADE_13[1]);
        multiplyFlat(a, tmp2, u, dimension);

        addScaledInPlaceFlat(v, a6, PADE_13[12]);
        addScaledInPlaceFlat(v, a4, PADE_13[10]);
        addScaledInPlaceFlat(v, a2, PADE_13[8]);
        multiplyFlat(a6, v, uInner, dimension);
        Arrays.fill(v, 0.0);
        addScaledInPlaceFlat(v, uInner, 1.0);
        addScaledInPlaceFlat(v, a6, PADE_13[6]);
        addScaledInPlaceFlat(v, a4, PADE_13[4]);
        addScaledInPlaceFlat(v, a2, PADE_13[2]);
        addScaledInPlaceFlat(v, identity, PADE_13[0]);

        final double[] numerator = new double[length];
        final double[] denominator = new double[length];
        addFlat(v, u, numerator, length);
        subtractFlat(v, u, denominator, length);
        solveFlat(denominator, 0, numerator, 0, out, outOffset, dimension, dimension);

        final double[] squared = new double[length];
        for (int s = 0; s < squarings; ++s) {
            multiplyFlat(out, outOffset, out, outOffset, squared, 0, dimension);
            System.arraycopy(squared, 0, out, outOffset, length);
        }
    }

    /**
     * Computes the Frechet derivative of the matrix exponential via the block identity
     * <pre>
     *   exp([[X, E], [0, X]]) = [[exp(X), L_exp(X, E)], [0, exp(X)]].
     * </pre>
     */
    public static void frechetExpFlat(final double[] x,
                                      final double[] direction,
                                      final double[] out,
                                      final int dimension) {
        frechetExpFlat(x, 0, direction, 0, out, 0, dimension);
    }

    public static void frechetExpFlat(final double[] x,
                                      final int xOffset,
                                      final double[] direction,
                                      final int directionOffset,
                                      final double[] out,
                                      final int outOffset,
                                      final int dimension) {
        final int twoN = 2 * dimension;
        final double[] block = new double[twoN * twoN];
        for (int i = 0; i < dimension; ++i) {
            final int xRow = xOffset + i * dimension;
            final int directionRow = directionOffset + i * dimension;
            final int topRow = i * twoN;
            final int bottomRow = (i + dimension) * twoN;
            for (int j = 0; j < dimension; ++j) {
                block[topRow + j] = x[xRow + j];
                block[topRow + j + dimension] = direction[directionRow + j];
                block[bottomRow + j + dimension] = x[xRow + j];
            }
        }

        final double[] expBlock = new double[twoN * twoN];
        expmFlat(block, expBlock, twoN);
        for (int i = 0; i < dimension; ++i) {
            final int outRow = outOffset + i * dimension;
            final int expRow = i * twoN;
            for (int j = 0; j < dimension; ++j) {
                out[outRow + j] = expBlock[expRow + j + dimension];
            }
        }
    }

    /**
     * Backpropagates through {@code Y = expm(X)}.
     *
     * <p>Given the upstream sensitivity {@code upstream = dL/dY}, returns
     * {@code out = dL/dX}. This is the Frobenius adjoint of the Frechet derivative,
     * computed as
     * <pre>
     *   dL/dX = L_exp(X^T, upstream)^T.
     * </pre>
     */
    public static void adjointExpFlat(final double[] x,
                                      final double[] upstream,
                                      final double[] out,
                                      final int dimension) {
        adjointExpFlat(x, 0, upstream, 0, out, 0, dimension);
    }

    public static void adjointExpFlat(final double[] x,
                                      final int xOffset,
                                      final double[] upstream,
                                      final int upstreamOffset,
                                      final double[] out,
                                      final int outOffset,
                                      final int dimension) {
        final int length = dimension * dimension;
        final double[] xTranspose = new double[length];
        final double[] tmp = new double[length];
        transposeFlat(x, xOffset, xTranspose, 0, dimension);
        frechetExpFlat(xTranspose, 0, upstream, upstreamOffset, tmp, 0, dimension);
        transposeFlat(tmp, 0, out, outOffset, dimension);
    }

    private static void addFlat(final double[] left,
                                final double[] right,
                                final double[] out,
                                final int length) {
        for (int i = 0; i < length; ++i) {
            out[i] = left[i] + right[i];
        }
    }

    private static void subtractFlat(final double[] left,
                                     final double[] right,
                                     final double[] out,
                                     final int length) {
        for (int i = 0; i < length; ++i) {
            out[i] = left[i] - right[i];
        }
    }

    private static void addScaledInPlaceFlat(final double[] accumulator,
                                             final double[] increment,
                                             final double scale) {
        for (int i = 0; i < accumulator.length; ++i) {
            accumulator[i] += scale * increment[i];
        }
    }

    private static void solveFlat(final double[] left,
                                  final int leftOffset,
                                  final double[] right,
                                  final int rightOffset,
                                  final double[] out,
                                  final int outOffset,
                                  final int n,
                                  final int m) {
        final double[] a = new double[n * n];
        final double[] b = new double[n * m];
        System.arraycopy(left, leftOffset, a, 0, n * n);
        System.arraycopy(right, rightOffset, b, 0, n * m);

        for (int pivot = 0; pivot < n; ++pivot) {
            int best = pivot;
            double bestAbs = Math.abs(a[pivot * n + pivot]);
            for (int row = pivot + 1; row < n; ++row) {
                final double candidate = Math.abs(a[row * n + pivot]);
                if (candidate > bestAbs) {
                    best = row;
                    bestAbs = candidate;
                }
            }
            if (!(bestAbs > 0.0)) {
                throw new IllegalArgumentException("Matrix is singular to working precision");
            }
            if (best != pivot) {
                swapRows(a, n, pivot, best);
                swapRows(b, m, pivot, best);
            }

            final int pivotRowA = pivot * n;
            final int pivotRowB = pivot * m;
            final double diagonal = a[pivotRowA + pivot];
            for (int j = pivot; j < n; ++j) {
                a[pivotRowA + j] /= diagonal;
            }
            for (int j = 0; j < m; ++j) {
                b[pivotRowB + j] /= diagonal;
            }

            for (int row = 0; row < n; ++row) {
                if (row == pivot) {
                    continue;
                }
                final int rowA = row * n;
                final int rowB = row * m;
                final double factor = a[rowA + pivot];
                if (factor == 0.0) {
                    continue;
                }
                for (int j = pivot; j < n; ++j) {
                    a[rowA + j] -= factor * a[pivotRowA + j];
                }
                for (int j = 0; j < m; ++j) {
                    b[rowB + j] -= factor * b[pivotRowB + j];
                }
            }
        }

        System.arraycopy(b, 0, out, outOffset, n * m);
    }

    private static void swapRows(final double[] matrix,
                                 final int columns,
                                 final int rowA,
                                 final int rowB) {
        final int offsetA = rowA * columns;
        final int offsetB = rowB * columns;
        for (int j = 0; j < columns; ++j) {
            final double tmp = matrix[offsetA + j];
            matrix[offsetA + j] = matrix[offsetB + j];
            matrix[offsetB + j] = tmp;
        }
    }

    private static double log2(final double x) {
        return Math.log(x) / Math.log(2.0);
    }
}
