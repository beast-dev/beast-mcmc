package dr.evomodel.continuous.ou;

/**
 * Dense matrix helpers for exact OU transitions.
 *
 * <p>This utility is intentionally self-contained so exact OU transition logic can live
 * without adding a new linear-algebra dependency.
 * The matrix exponential uses scaling-and-squaring with the order-13 Pad\'e approximant.
 *
 * <p>For selection-matrix gradients we expose an <em>adjoint</em> (backpropagation) helper
 * for the matrix exponential.  If
 * <pre>
 *     Y = expm(X)
 * </pre>
 * and the upstream sensitivity is G = dL/dY, then the Frobenius-adjoint sensitivity wrt X is
 * <pre>
 *     dL/dX = L_exp(X^T, G)^T,
 * </pre>
 * where L_exp is the Fr\'echet derivative of the matrix exponential.  The implementation
 * evaluates this adjoint via the standard block-exponential identity rather than by explicit
 * forward-direction loops over basis matrices.
 */
public final class MatrixExponentialUtils {

    private static final double THETA_13 = 5.371920351148152;

    private MatrixExponentialUtils() {
        // no instances
    }

    public static void setIdentity(final double[][] matrix) {
        for (int i = 0; i < matrix.length; ++i) {
            for (int j = 0; j < matrix[i].length; ++j) {
                matrix[i][j] = (i == j) ? 1.0 : 0.0;
            }
        }
    }

    public static void fill(final double[][] matrix, final double value) {
        for (int i = 0; i < matrix.length; ++i) {
            for (int j = 0; j < matrix[i].length; ++j) {
                matrix[i][j] = value;
            }
        }
    }

    public static void copy(final double[][] source, final double[][] target) {
        for (int i = 0; i < source.length; ++i) {
            System.arraycopy(source[i], 0, target[i], 0, source[i].length);
        }
    }

    public static void transpose(final double[][] source, final double[][] out) {
        for (int i = 0; i < source.length; ++i) {
            for (int j = 0; j < source[i].length; ++j) {
                out[j][i] = source[i][j];
            }
        }
    }

    public static void scale(final double[][] matrix, final double scalar, final double[][] out) {
        for (int i = 0; i < matrix.length; ++i) {
            for (int j = 0; j < matrix[i].length; ++j) {
                out[i][j] = scalar * matrix[i][j];
            }
        }
    }

    public static void add(final double[][] left, final double[][] right, final double[][] out) {
        for (int i = 0; i < left.length; ++i) {
            for (int j = 0; j < left[i].length; ++j) {
                out[i][j] = left[i][j] + right[i][j];
            }
        }
    }

    public static void subtract(final double[][] left, final double[][] right, final double[][] out) {
        for (int i = 0; i < left.length; ++i) {
            for (int j = 0; j < left[i].length; ++j) {
                out[i][j] = left[i][j] - right[i][j];
            }
        }
    }

    public static void addScaledInPlace(final double[][] accumulator, final double[][] increment, final double scale) {
        for (int i = 0; i < accumulator.length; ++i) {
            for (int j = 0; j < accumulator[i].length; ++j) {
                accumulator[i][j] += scale * increment[i][j];
            }
        }
    }

    public static void multiply(final double[][] left, final double[][] right, final double[][] out) {
        final int rows = left.length;
        final int inner = right.length;
        final int cols = right[0].length;
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                double sum = 0.0;
                for (int k = 0; k < inner; ++k) {
                    sum += left[i][k] * right[k][j];
                }
                out[i][j] = sum;
            }
        }
    }

    public static void multiply(final double[][] matrix, final double[] vector, final double[] out) {
        for (int i = 0; i < matrix.length; ++i) {
            double sum = 0.0;
            for (int j = 0; j < vector.length; ++j) {
                sum += matrix[i][j] * vector[j];
            }
            out[i] = sum;
        }
    }

    public static void outerProduct(final double[] left, final double[] right, final double[][] out) {
        for (int i = 0; i < left.length; ++i) {
            for (int j = 0; j < right.length; ++j) {
                out[i][j] = left[i] * right[j];
            }
        }
    }

    public static void symmetrize(final double[][] matrix) {
        for (int i = 0; i < matrix.length; ++i) {
            for (int j = i + 1; j < matrix[i].length; ++j) {
                final double average = 0.5 * (matrix[i][j] + matrix[j][i]);
                matrix[i][j] = average;
                matrix[j][i] = average;
            }
        }
    }

    public static void transposeFlat(final double[] source, final double[] out, final int dimension) {
        for (int i = 0; i < dimension; ++i) {
            final int rowOffset = i * dimension;
            for (int j = 0; j < dimension; ++j) {
                out[j * dimension + i] = source[rowOffset + j];
            }
        }
    }

    public static void multiplyFlat(final double[] left,
                                    final double[] right,
                                    final double[] out,
                                    final int dimension) {
        for (int i = 0; i < dimension; ++i) {
            final int rowOffset = i * dimension;
            for (int j = 0; j < dimension; ++j) {
                double sum = 0.0;
                for (int k = 0; k < dimension; ++k) {
                    sum += left[rowOffset + k] * right[k * dimension + j];
                }
                out[rowOffset + j] = sum;
            }
        }
    }

    public static void symmetrizeFlat(final double[] matrix, final int dimension) {
        for (int i = 0; i < dimension; ++i) {
            for (int j = i + 1; j < dimension; ++j) {
                final int ij = i * dimension + j;
                final int ji = j * dimension + i;
                final double average = 0.5 * (matrix[ij] + matrix[ji]);
                matrix[ij] = average;
                matrix[ji] = average;
            }
        }
    }

    public static void expmFlat(final double[] matrix, final double[] out, final int dimension) {
        final double[][] dense = new double[dimension][dimension];
        final double[][] denseOut = new double[dimension][dimension];
        for (int i = 0; i < dimension; ++i) {
            System.arraycopy(matrix, i * dimension, dense[i], 0, dimension);
        }
        expm(dense, denseOut);
        for (int i = 0; i < dimension; ++i) {
            System.arraycopy(denseOut[i], 0, out, i * dimension, dimension);
        }
    }

    public static double norm1(final double[][] matrix) {
        double best = 0.0;
        for (int j = 0; j < matrix[0].length; ++j) {
            double sum = 0.0;
            for (double[] doubles : matrix) {
                sum += Math.abs(doubles[j]);
            }
            if (sum > best) {
                best = sum;
            }
        }
        return best;
    }

    public static void solve(final double[][] left, final double[][] right, final double[][] out) {
        final int n = left.length;
        final int m = right[0].length;
        final double[][] a = new double[n][n];
        final double[][] b = new double[n][m];
        copy(left, a);
        copy(right, b);

        for (int pivot = 0; pivot < n; ++pivot) {
            int best = pivot;
            double bestAbs = Math.abs(a[pivot][pivot]);
            for (int row = pivot + 1; row < n; ++row) {
                final double candidate = Math.abs(a[row][pivot]);
                if (candidate > bestAbs) {
                    best = row;
                    bestAbs = candidate;
                }
            }
            if (!(bestAbs > 0.0)) {
                throw new IllegalArgumentException("Matrix is singular to working precision");
            }
            if (best != pivot) {
                final double[] tmpA = a[pivot];
                a[pivot] = a[best];
                a[best] = tmpA;
                final double[] tmpB = b[pivot];
                b[pivot] = b[best];
                b[best] = tmpB;
            }

            final double diagonal = a[pivot][pivot];
            for (int j = pivot; j < n; ++j) {
                a[pivot][j] /= diagonal;
            }
            for (int j = 0; j < m; ++j) {
                b[pivot][j] /= diagonal;
            }

            for (int row = 0; row < n; ++row) {
                if (row == pivot) {
                    continue;
                }
                final double factor = a[row][pivot];
                if (factor == 0.0) {
                    continue;
                }
                for (int j = pivot; j < n; ++j) {
                    a[row][j] -= factor * a[pivot][j];
                }
                for (int j = 0; j < m; ++j) {
                    b[row][j] -= factor * b[pivot][j];
                }
            }
        }

        copy(b, out);
    }

    public static void expm(final double[][] matrix, final double[][] out) {
        final int n = matrix.length;
        final double[][] a = new double[n][n];
        copy(matrix, a);

        final double aNorm = norm1(a);
        final int squarings = aNorm <= THETA_13 ? 0 : Math.max(0, (int) Math.ceil(log2(aNorm / THETA_13)));
        if (squarings > 0) {
            final double scale = Math.scalb(1.0, -squarings);
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    a[i][j] *= scale;
                }
            }
        }

        final double[][] a2 = new double[n][n];
        final double[][] a4 = new double[n][n];
        final double[][] a6 = new double[n][n];
        multiply(a, a, a2);
        multiply(a2, a2, a4);
        multiply(a2, a4, a6);

        final double[] b = {
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

        final double[][] tmp1 = new double[n][n];
        final double[][] tmp2 = new double[n][n];
        final double[][] uInner = new double[n][n];
        final double[][] u = new double[n][n];
        final double[][] v = new double[n][n];
        final double[][] identity = new double[n][n];
        setIdentity(identity);

        fill(tmp1, 0.0);
        addScaledInPlace(tmp1, a6, b[13]);
        addScaledInPlace(tmp1, a4, b[11]);
        addScaledInPlace(tmp1, a2, b[9]);
        multiply(a6, tmp1, tmp2);
        addScaledInPlace(tmp2, a6, b[7]);
        addScaledInPlace(tmp2, a4, b[5]);
        addScaledInPlace(tmp2, a2, b[3]);
        addScaledInPlace(tmp2, identity, b[1]);
        multiply(a, tmp2, u);

        fill(v, 0.0);
        addScaledInPlace(v, a6, b[12]);
        addScaledInPlace(v, a4, b[10]);
        addScaledInPlace(v, a2, b[8]);
        multiply(a6, v, uInner);
        fill(v, 0.0);
        addScaledInPlace(v, uInner, 1.0);
        addScaledInPlace(v, a6, b[6]);
        addScaledInPlace(v, a4, b[4]);
        addScaledInPlace(v, a2, b[2]);
        addScaledInPlace(v, identity, b[0]);

        final double[][] numerator = new double[n][n];
        final double[][] denominator = new double[n][n];
        add(v, u, numerator);
        subtract(v, u, denominator);
        solve(denominator, numerator, out);

        final double[][] squared = new double[n][n];
        for (int s = 0; s < squarings; ++s) {
            multiply(out, out, squared);
            copy(squared, out);
        }
    }

    /**
     * Computes the Fr\'echet derivative of the matrix exponential via the block identity
     * <pre>
     *   exp([[X, E], [0, X]]) = [[exp(X), L_exp(X, E)], [0, exp(X)]].
     * </pre>
     */
    public static void frechetExp(final double[][] x, final double[][] direction, final double[][] out) {
        final int n = x.length;
        final int twoN = 2 * n;
        final double[][] block = new double[twoN][twoN];
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                block[i][j] = x[i][j];
                block[i][j + n] = direction[i][j];
                block[i + n][j] = 0.0;
                block[i + n][j + n] = x[i][j];
            }
        }
        final double[][] expBlock = new double[twoN][twoN];
        expm(block, expBlock);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                out[i][j] = expBlock[i][j + n];
            }
        }
    }

    /**
     * Backpropagates through {@code Y = expm(X)}.
     *
     * <p>Given the upstream sensitivity {@code upstream = dL/dY}, returns
     * {@code out = dL/dX}.  This is the Frobenius adjoint of the Fr\'echet derivative,
     * computed as
     * <pre>
     *   dL/dX = L_exp(X^T, upstream)^T.
     * </pre>
     */
    public static void adjointExp(final double[][] x, final double[][] upstream, final double[][] out) {
        final int n = x.length;
        final double[][] xTranspose = new double[n][n];
        final double[][] tmp = new double[n][n];
        transpose(x, xTranspose);
        frechetExp(xTranspose, upstream, tmp);
        transpose(tmp, out);
    }

    public static void adjointExpFlat(final double[] x,
                                      final double[] upstream,
                                      final double[] out,
                                      final int dimension) {
        final double[][] denseX = new double[dimension][dimension];
        final double[][] denseUpstream = new double[dimension][dimension];
        final double[][] denseOut = new double[dimension][dimension];
        for (int i = 0; i < dimension; ++i) {
            System.arraycopy(x, i * dimension, denseX[i], 0, dimension);
            System.arraycopy(upstream, i * dimension, denseUpstream[i], 0, dimension);
        }
        adjointExp(denseX, denseUpstream, denseOut);
        for (int i = 0; i < dimension; ++i) {
            System.arraycopy(denseOut[i], 0, out, i * dimension, dimension);
        }
    }

    private static double log2(final double x) {
        return Math.log(x) / Math.log(2.0);
    }
}
