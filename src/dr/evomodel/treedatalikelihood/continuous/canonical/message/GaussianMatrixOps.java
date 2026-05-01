package dr.evomodel.treedatalikelihood.continuous.canonical.message;

/**
 * Shared dense linear-algebra utilities for Gaussian likelihood engines.
 *
 * <p>This helper centralizes the matrix and Cholesky routines that are shared by
 * expectation-form and canonical-form Gaussian engines.
 */
public final class GaussianMatrixOps {

    public static final double LOG_TWO_PI = Math.log(2.0 * Math.PI);
    public static final double MIN_DIAGONAL_JITTER = 1E-10;

    private static final double SYMMETRY_TOLERANCE = 1E-12;

    private GaussianMatrixOps() { }

    public static void multiplyMatrixVector(final double[][] matrix,
                                     final double[] vector,
                                     final double[] out) {
        multiplyMatrixVector(matrix, vector, out, matrix.length, vector.length);
    }

    public static void multiplyMatrixVector(final double[][] matrix,
                                     final double[] vector,
                                     final double[] out,
                                     final int rows,
                                     final int cols) {
        for (int i = 0; i < rows; ++i) {
            double sum = 0.0;
            for (int j = 0; j < cols; ++j) {
                sum += matrix[i][j] * vector[j];
            }
            out[i] = sum;
        }
    }

    public static void multiplyMatrixMatrix(final double[][] left,
                                     final double[][] right,
                                     final double[][] out) {
        multiplyMatrixMatrix(left, right, out, left.length, right.length, right[0].length);
    }

    public static void multiplyMatrixMatrix(final double[][] left,
                                     final double[][] right,
                                     final double[][] out,
                                     final int rows,
                                     final int inner,
                                     final int cols) {
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

    public static void multiplyMatrixMatrixTransposedRight(final double[][] left,
                                                    final double[][] right,
                                                    final double[][] out) {
        final int rows = left.length;
        final int inner = left[0].length;
        final int cols = right.length;
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                double sum = 0.0;
                for (int k = 0; k < inner; ++k) {
                    sum += left[i][k] * right[j][k];
                }
                out[i][j] = sum;
            }
        }
    }

    public static void addMatrixInPlace(final double[][] accumulator, final double[][] increment) {
        for (int i = 0; i < accumulator.length; ++i) {
            for (int j = 0; j < accumulator[i].length; ++j) {
                accumulator[i][j] += increment[i][j];
            }
        }
    }

    public static void addInPlace(final double[] accumulator, final double[] increment) {
        for (int i = 0; i < accumulator.length; ++i) {
            accumulator[i] += increment[i];
        }
    }

    // -----------------------------------------------------------------------
    // Flat (row-major) variants – index convention: matrix[i * dim + j]
    // -----------------------------------------------------------------------

    public static void multiplyMatrixMatrixFlat(final double[] left,
                                         final double[] right,
                                         final double[] out,
                                         final int dim) {
        for (int i = 0; i < dim; ++i) {
            final int iOff = i * dim;
            for (int j = 0; j < dim; ++j) {
                double sum = 0.0;
                for (int k = 0; k < dim; ++k) {
                    sum += left[iOff + k] * right[k * dim + j];
                }
                out[iOff + j] = sum;
            }
        }
    }

    public static void multiplyMatrixVectorFlat(final double[] matrix,
                                         final double[] vector,
                                         final double[] out,
                                         final int dim) {
        for (int i = 0; i < dim; ++i) {
            final int iOff = i * dim;
            double sum = 0.0;
            for (int j = 0; j < dim; ++j) {
                sum += matrix[iOff + j] * vector[j];
            }
            out[i] = sum;
        }
    }

    public static void symmetrizeFlat(final double[] matrix, final int dim) {
        for (int i = 0; i < dim; ++i) {
            for (int j = i + 1; j < dim; ++j) {
                final double avg = 0.5 * (matrix[i * dim + j] + matrix[j * dim + i]);
                matrix[i * dim + j] = avg;
                matrix[j * dim + i] = avg;
            }
        }
    }

    public static void copyMatrixFlat(final double[] source, final double[] target, final int dim) {
        System.arraycopy(source, 0, target, 0, dim * dim);
    }

    /** Copy from a 2-D row-major {@code double[][]} into a flat row-major {@code double[]}. */
    public static void copyMatrixToFlat(final double[][] source, final double[] target, final int dim) {
        matrixToRowMajor(source, target, dim);
    }

    /** Copy from a 2-D row-major {@code double[][]} into a flat row-major {@code double[]}. */
    public static void matrixToRowMajor(final double[][] source, final double[] target, final int dim) {
        for (int i = 0; i < dim; ++i) {
            System.arraycopy(source[i], 0, target, i * dim, dim);
        }
    }

    /** Copy from a flat row-major {@code double[]} into a 2-D row-major {@code double[][]}. */
    public static void copyFlatToMatrix(final double[] source, final double[][] target, final int dim) {
        rowMajorToMatrix(source, target, dim);
    }

    /** Copy from a flat row-major {@code double[]} into a 2-D row-major {@code double[][]}. */
    public static void rowMajorToMatrix(final double[] source, final double[][] target, final int dim) {
        for (int i = 0; i < dim; ++i) {
            System.arraycopy(source, i * dim, target[i], 0, dim);
        }
    }

    /** Copy a flat row-major square matrix into flat column-major parameter order. */
    public static void rowMajorToColumnMajorParameter(final double[] source,
                                                      final double[] target,
                                                      final int dim) {
        for (int row = 0; row < dim; ++row) {
            for (int col = 0; col < dim; ++col) {
                target[col * dim + row] = source[row * dim + col];
            }
        }
    }

    /** Copy a 2-D row-major square matrix into flat column-major parameter order. */
    public static void matrixToColumnMajorParameter(final double[][] source,
                                                    final double[] target,
                                                    final int dim) {
        int index = 0;
        for (int col = 0; col < dim; ++col) {
            for (int row = 0; row < dim; ++row) {
                target[index++] = source[row][col];
            }
        }
    }

    /** Copy a flat column-major parameter vector into flat row-major square-matrix order. */
    public static void columnMajorParameterToRowMajor(final double[] source,
                                                      final double[] target,
                                                      final int dim) {
        for (int row = 0; row < dim; ++row) {
            for (int col = 0; col < dim; ++col) {
                target[row * dim + col] = source[col * dim + row];
            }
        }
    }

    /** Copy the transpose of a flat row-major square matrix into a 2-D matrix. */
    public static void transposeFlatToMatrix(final double[] source, final double[][] target, final int dim) {
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                target[j][i] = source[i * dim + j];
            }
        }
    }

    public static void transposeFlatSquareInPlace(final double[] matrix, final int dim) {
        for (int i = 0; i < dim; ++i) {
            for (int j = i + 1; j < dim; ++j) {
                final int ij = i * dim + j;
                final int ji = j * dim + i;
                final double tmp = matrix[ij];
                matrix[ij] = matrix[ji];
                matrix[ji] = tmp;
            }
        }
    }

    public static void multiplyFlatByMatrix(final double[] left, final double[][] right,
                                            final double[][] out, final int dim) {
        for (int i = 0; i < dim; ++i) {
            final int iOff = i * dim;
            for (int j = 0; j < dim; ++j) {
                double sum = 0.0;
                for (int k = 0; k < dim; ++k) {
                    sum += left[iOff + k] * right[k][j];
                }
                out[i][j] = sum;
            }
        }
    }

    public static void multiplyMatrixByFlat(final double[][] left, final double[] right,
                                            final double[][] out, final int dim) {
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                double sum = 0.0;
                for (int k = 0; k < dim; ++k) {
                    sum += left[i][k] * right[k * dim + j];
                }
                out[i][j] = sum;
            }
        }
    }

    /**
     * Cholesky factorization of a flat row-major symmetric positive-definite matrix.
     * Returns a {@link FlatCholeskyFactor}.
     */
    public static FlatCholeskyFactor choleskyFlat(final double[] matrix, final int dim) {
        final double[] lower = new double[dim * dim];
        if (!tryCholeskyFlat(matrix, lower, dim)) {
            throw new IllegalArgumentException("Matrix is not positive definite");
        }
        return new FlatCholeskyFactor(lower, dim);
    }

    public static boolean tryCholeskyFlat(final double[] matrix,
                                    final double[] lowerOut,
                                    final int dim) {
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j <= i; ++j) {
                double sum = matrix[i * dim + j];
                for (int k = 0; k < j; ++k) {
                    sum -= lowerOut[i * dim + k] * lowerOut[j * dim + k];
                }
                if (i == j) {
                    if (sum <= 0.0) return false;
                    lowerOut[i * dim + j] = Math.sqrt(sum);
                } else {
                    final double denom = lowerOut[j * dim + j];
                    if (denom == 0.0) return false;
                    lowerOut[i * dim + j] = sum / denom;
                }
            }
            for (int j = i + 1; j < dim; ++j) {
                lowerOut[i * dim + j] = 0.0;
            }
        }
        return true;
    }

    public static void invertPositiveDefiniteFromFlatCholesky(final double[] out,
                                                        final FlatCholeskyFactor factor) {
        final int dim = factor.dimension;
        final double[] lowerInverse = new double[dim * dim];
        final double[] solution = new double[dim];
        final double[] basis = new double[dim];

        for (int column = 0; column < dim; ++column) {
            for (int i = 0; i < dim; ++i) basis[i] = 0.0;
            basis[column] = 1.0;
            factor.solveSymmetricSystem(basis, solution);
            for (int row = 0; row < dim; ++row) {
                lowerInverse[row * dim + column] = solution[row];
            }
        }
        System.arraycopy(lowerInverse, 0, out, 0, dim * dim);
        symmetrizeFlat(out, dim);
    }

    public static final class FlatCholeskyFactor {
        final double[] lower;
        final int dimension;

        FlatCholeskyFactor(final double[] lower, final int dimension) {
            this.lower = lower;
            this.dimension = dimension;
        }

        public double logDeterminant() {
            double value = 0.0;
            for (int i = 0; i < dimension; ++i) {
                value += 2.0 * Math.log(lower[i * dimension + i]);
            }
            return value;
        }

        public void solveSymmetricSystem(final double[] rhs, final double[] out) {
            final int dim = dimension;
            final double[] y = new double[dim];
            for (int i = 0; i < dim; ++i) {
                double sum = rhs[i];
                for (int j = 0; j < i; ++j) {
                    sum -= lower[i * dim + j] * y[j];
                }
                y[i] = sum / lower[i * dim + i];
            }
            for (int i = dim - 1; i >= 0; --i) {
                double sum = y[i];
                for (int j = i + 1; j < dim; ++j) {
                    sum -= lower[j * dim + i] * out[j];
                }
                out[i] = sum / lower[i * dim + i];
            }
        }
    }

    public static void identityMinus(final double[][] matrix) {
        for (int i = 0; i < matrix.length; ++i) {
            for (int j = 0; j < matrix[i].length; ++j) {
                matrix[i][j] = (i == j ? 1.0 : 0.0) - matrix[i][j];
            }
        }
    }

    public static void copyVector(final double[] source, final double[] target) {
        System.arraycopy(source, 0, target, 0, source.length);
    }

    public static void copyMatrix(final double[][] source, final double[][] target) {
        copyMatrix(source, target, source.length, source[0].length);
    }

    public static void copyMatrix(final double[][] source,
                           final double[][] target,
                           final int rows,
                           final int cols) {
        for (int i = 0; i < rows; ++i) {
            System.arraycopy(source[i], 0, target[i], 0, cols);
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

    public static void addJitterToDiagonal(final double[][] matrix, final double minimumJitter) {
        for (int i = 0; i < matrix.length; ++i) {
            if (matrix[i][i] < minimumJitter) {
                matrix[i][i] = minimumJitter;
            }
        }
    }

    public static double quadraticForm(final double[][] matrix, final double[] vector) {
        double result = 0.0;
        for (int i = 0; i < vector.length; ++i) {
            double rowSum = 0.0;
            for (int j = 0; j < vector.length; ++j) {
                rowSum += matrix[i][j] * vector[j];
            }
            result += vector[i] * rowSum;
        }
        return result;
    }

    public static CholeskyFactor cholesky(final double[][] matrix) {
        final int dim = matrix.length;
        final double[][] lower = new double[dim][dim];

        if (!tryCholesky(matrix, lower, dim)) {
            throw new IllegalArgumentException("Matrix is not positive definite");
        }

        return new CholeskyFactor(lower);
    }

    public static boolean tryCholesky(final double[][] matrix,
                               final double[][] lowerOut,
                               final int dimension) {
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j <= i; ++j) {
                double sum = matrix[i][j];
                for (int k = 0; k < j; ++k) {
                    sum -= lowerOut[i][k] * lowerOut[j][k];
                }

                if (i == j) {
                    if (sum <= 0.0) {
                        return false;
                    }
                    lowerOut[i][j] = Math.sqrt(sum);
                } else {
                    final double denom = lowerOut[j][j];
                    if (denom == 0.0) {
                        return false;
                    }
                    lowerOut[i][j] = sum / denom;
                }
            }

            for (int j = i + 1; j < dimension; ++j) {
                if (Math.abs(matrix[i][j] - matrix[j][i]) > SYMMETRY_TOLERANCE) {
                    return false;
                }
            }
            for (int j = i + 1; j < dimension; ++j) {
                lowerOut[i][j] = 0.0;
            }
        }
        return true;
    }

    public static void invertPositiveDefiniteFromCholesky(final double[][] out,
                                                   final CholeskyFactor factor) {
        final int dim = out.length;
        final double[][] inverse = new double[dim][dim];
        final double[] basis = new double[dim];
        final double[] solution = new double[dim];

        for (int column = 0; column < dim; ++column) {
            for (int i = 0; i < dim; ++i) {
                basis[i] = 0.0;
            }
            basis[column] = 1.0;
            factor.solveSymmetricSystem(basis, solution);
            for (int row = 0; row < dim; ++row) {
                inverse[row][column] = solution[row];
            }
        }

        copyMatrix(inverse, out);
        symmetrize(out);
    }

    public static void invertPositiveDefiniteFromLowerTriangular(final double[][] inverseOut,
                                                          final double[][] lowerTriangular,
                                                          final double[] yScratch,
                                                          final int dimension) {
        for (int column = 0; column < dimension; ++column) {
            for (int i = 0; i < dimension; ++i) {
                double sum = (i == column) ? 1.0 : 0.0;
                for (int j = 0; j < i; ++j) {
                    sum -= lowerTriangular[i][j] * yScratch[j];
                }
                yScratch[i] = sum / lowerTriangular[i][i];
            }

            for (int i = dimension - 1; i >= 0; --i) {
                double sum = yScratch[i];
                for (int j = i + 1; j < dimension; ++j) {
                    sum -= lowerTriangular[j][i] * inverseOut[j][column];
                }
                inverseOut[i][column] = sum / lowerTriangular[i][i];
            }
        }
        symmetrize(inverseOut);
    }

    public static final class CholeskyFactor {
        private final double[][] lower;

        CholeskyFactor(final double[][] lower) {
            this.lower = lower;
        }

        public double logDeterminant() {
            double value = 0.0;
            for (int i = 0; i < lower.length; ++i) {
                value += 2.0 * Math.log(lower[i][i]);
            }
            return value;
        }

        public void solveSymmetricSystem(final double[] rhs, final double[] out) {
            final int dim = lower.length;
            final double[] y = new double[dim];

            for (int i = 0; i < dim; ++i) {
                double sum = rhs[i];
                for (int j = 0; j < i; ++j) {
                    sum -= lower[i][j] * y[j];
                }
                y[i] = sum / lower[i][i];
            }

            for (int i = dim - 1; i >= 0; --i) {
                double sum = y[i];
                for (int j = i + 1; j < dim; ++j) {
                    sum -= lower[j][i] * out[j];
                }
                out[i] = sum / lower[i][i];
            }
        }
    }
}
