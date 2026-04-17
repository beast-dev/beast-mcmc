package dr.inference.timeseries.engine.gaussian;

/**
 * Shared dense linear-algebra utilities for Gaussian likelihood engines.
 *
 * <p>This helper centralizes the matrix and Cholesky routines that are shared by
 * expectation-form and canonical-form Gaussian engines. Keeping them here avoids
 * coupling canonical code to {@link KalmanLikelihoodEngine} by name.
 */
final class GaussianMatrixOps {

    static final double LOG_TWO_PI = Math.log(2.0 * Math.PI);
    static final double MIN_DIAGONAL_JITTER = 1E-10;

    private static final double SYMMETRY_TOLERANCE = 1E-12;

    private GaussianMatrixOps() { }

    static void multiplyMatrixVector(final double[][] matrix,
                                     final double[] vector,
                                     final double[] out) {
        multiplyMatrixVector(matrix, vector, out, matrix.length, vector.length);
    }

    static void multiplyMatrixVector(final double[][] matrix,
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

    static void multiplyMatrixMatrix(final double[][] left,
                                     final double[][] right,
                                     final double[][] out) {
        multiplyMatrixMatrix(left, right, out, left.length, right.length, right[0].length);
    }

    static void multiplyMatrixMatrix(final double[][] left,
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

    static void multiplyMatrixMatrixTransposedRight(final double[][] left,
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

    static void addMatrixInPlace(final double[][] accumulator, final double[][] increment) {
        for (int i = 0; i < accumulator.length; ++i) {
            for (int j = 0; j < accumulator[i].length; ++j) {
                accumulator[i][j] += increment[i][j];
            }
        }
    }

    static void addInPlace(final double[] accumulator, final double[] increment) {
        for (int i = 0; i < accumulator.length; ++i) {
            accumulator[i] += increment[i];
        }
    }

    static void identityMinus(final double[][] matrix) {
        for (int i = 0; i < matrix.length; ++i) {
            for (int j = 0; j < matrix[i].length; ++j) {
                matrix[i][j] = (i == j ? 1.0 : 0.0) - matrix[i][j];
            }
        }
    }

    static void copyVector(final double[] source, final double[] target) {
        System.arraycopy(source, 0, target, 0, source.length);
    }

    static void copyMatrix(final double[][] source, final double[][] target) {
        copyMatrix(source, target, source.length, source[0].length);
    }

    static void copyMatrix(final double[][] source,
                           final double[][] target,
                           final int rows,
                           final int cols) {
        for (int i = 0; i < rows; ++i) {
            System.arraycopy(source[i], 0, target[i], 0, cols);
        }
    }

    static void symmetrize(final double[][] matrix) {
        for (int i = 0; i < matrix.length; ++i) {
            for (int j = i + 1; j < matrix[i].length; ++j) {
                final double average = 0.5 * (matrix[i][j] + matrix[j][i]);
                matrix[i][j] = average;
                matrix[j][i] = average;
            }
        }
    }

    static void addJitterToDiagonal(final double[][] matrix, final double minimumJitter) {
        for (int i = 0; i < matrix.length; ++i) {
            if (matrix[i][i] < minimumJitter) {
                matrix[i][i] = minimumJitter;
            }
        }
    }

    static double quadraticForm(final double[][] matrix, final double[] vector) {
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

    static CholeskyFactor cholesky(final double[][] matrix) {
        final int dim = matrix.length;
        final double[][] lower = new double[dim][dim];

        if (!tryCholesky(matrix, lower, dim)) {
            throw new IllegalArgumentException("Matrix is not positive definite");
        }

        return new CholeskyFactor(lower);
    }

    static boolean tryCholesky(final double[][] matrix,
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

    static void invertPositiveDefiniteFromCholesky(final double[][] out,
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

    static void invertPositiveDefiniteFromLowerTriangular(final double[][] inverseOut,
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

    static final class CholeskyFactor {
        private final double[][] lower;

        CholeskyFactor(final double[][] lower) {
            this.lower = lower;
        }

        double logDeterminant() {
            double value = 0.0;
            for (int i = 0; i < lower.length; ++i) {
                value += 2.0 * Math.log(lower[i][i]);
            }
            return value;
        }

        void solveSymmetricSystem(final double[] rhs, final double[] out) {
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
