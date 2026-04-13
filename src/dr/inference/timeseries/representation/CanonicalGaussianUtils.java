package dr.inference.timeseries.representation;

/**
 * Utilities for converting expectation-form Gaussian quantities into canonical form.
 */
public final class CanonicalGaussianUtils {

    private static final double SYMMETRIC_JITTER_RELATIVE = 1.0e-12;
    private static final double SYMMETRIC_JITTER_ABSOLUTE = 1.0e-12;
    private static final ThreadLocal<Workspace> WORKSPACE =
            ThreadLocal.withInitial(Workspace::new);

    private CanonicalGaussianUtils() {
        // no instances
    }

    public static void fillStateFromMoments(final double[] mean,
                                            final double[][] covariance,
                                            final CanonicalGaussianState out) {
        final int dimension = mean.length;
        validateSquare(covariance, dimension, "covariance");
        if (out.getDimension() != dimension) {
            throw new IllegalArgumentException("State output dimension mismatch");
        }

        final double[][] precision = out.precision;
        final Workspace workspace = workspaceFor(dimension);
        final double logDet = invertSymmetricPositiveDefinite(covariance, precision, dimension, workspace);
        multiply(precision, mean, out.information);
        out.logNormalizer =
                0.5 * (dimension * Math.log(2.0 * Math.PI) + logDet + dot(mean, out.information));
    }

    public static void fillTransitionFromMoments(final double[][] transitionMatrix,
                                                 final double[] transitionOffset,
                                                 final double[][] transitionCovariance,
                                                 final CanonicalGaussianTransition out) {
        final int dimension = transitionOffset.length;
        validateSquare(transitionMatrix, dimension, "transitionMatrix");
        validateSquare(transitionCovariance, dimension, "transitionCovariance");
        if (out.getDimension() != dimension) {
            throw new IllegalArgumentException("Transition output dimension mismatch");
        }

        final double[][] precision = new double[dimension][dimension];
        final Workspace workspace = workspaceFor(dimension);
        final double logDet = invertSymmetricPositiveDefinite(transitionCovariance, precision, dimension, workspace);
        final double[][] transitionTranspose = new double[dimension][dimension];
        final double[][] tmp = new double[dimension][dimension];
        transpose(transitionMatrix, transitionTranspose);

        multiply(transitionTranspose, precision, tmp);
        multiply(tmp, transitionMatrix, out.precisionXX);

        multiply(transitionTranspose, precision, out.precisionXY);
        scaleInPlace(out.precisionXY, -1.0);

        multiply(precision, transitionMatrix, out.precisionYX);
        scaleInPlace(out.precisionYX, -1.0);

        copy(precision, out.precisionYY);

        multiply(precision, transitionOffset, out.informationY);
        multiply(transitionTranspose, out.informationY, out.informationX);
        scaleInPlace(out.informationX, -1.0);

        out.logNormalizer =
                0.5 * (dimension * Math.log(2.0 * Math.PI)
                        + logDet
                        + dot(transitionOffset, out.informationY));
    }

    public static void fillMomentsFromCanonical(final CanonicalGaussianState state,
                                                final double[] meanOut,
                                                final double[][] covarianceOut) {
        final int dimension = state.getDimension();
        if (meanOut.length != dimension) {
            throw new IllegalArgumentException("Mean output dimension mismatch");
        }
        validateSquare(covarianceOut, dimension, "covarianceOut");

        final Workspace workspace = workspaceFor(dimension);
        try {
            // Fast path: avoid symmetrization/jitter unless needed.
            invertSymmetricPositiveDefinite(state.precision, covarianceOut, dimension, workspace);
        } catch (IllegalArgumentException ignored) {
            invertSymmetricPositiveDefiniteWithJitter(state.precision, covarianceOut, dimension, workspace);
        }
        multiply(covarianceOut, state.information, meanOut);
    }

    private static void invertSymmetricPositiveDefiniteWithJitter(final double[][] matrix,
                                                                  final double[][] inverseOut,
                                                                  final int dimension,
                                                                  final Workspace workspace) {
        final double[][] symmetric = workspace.symmetric;
        final double[][] adjusted = workspace.adjusted;
        symmetrizeCopy(matrix, symmetric, dimension);

        final double jitterBase = Math.max(
                SYMMETRIC_JITTER_ABSOLUTE,
                SYMMETRIC_JITTER_RELATIVE * Math.max(1.0, maxAbsDiagonal(symmetric, dimension)));

        double jitter = 0.0;
        double lowerBound = Double.NaN;
        for (int attempt = 0; attempt < 12; ++attempt) {
            copySquare(symmetric, adjusted, dimension);
            if (jitter > 0.0) {
                addDiagonalJitter(adjusted, jitter, dimension);
            }
            try {
                invertSymmetricPositiveDefinite(adjusted, inverseOut, dimension, workspace);
                return;
            } catch (IllegalArgumentException ignored) {
                // retry with larger diagonal jitter
            }
            if (jitter == 0.0) {
                if (Double.isNaN(lowerBound)) {
                    lowerBound = gershgorinLowerBound(symmetric, dimension);
                }
                jitter = lowerBound > 0.0 ? jitterBase : (-lowerBound + jitterBase);
            } else {
                jitter *= 10.0;
            }
        }

        throw new IllegalArgumentException(
                "Matrix is not positive definite (failed robust inversion retries)");
    }

    private static double invertSymmetricPositiveDefinite(final double[][] matrix,
                                                          final double[][] inverseOut,
                                                          final int dimension,
                                                          final Workspace workspace) {
        final double[][] cholesky = workspace.cholesky;
        final double[][] lowerInverse = workspace.lowerInverse;
        final double[][] lowerInverseTranspose = workspace.lowerInverseTranspose;
        copySquare(matrix, cholesky, dimension);
        double logDet = 0.0;

        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j <= i; ++j) {
                double sum = cholesky[i][j];
                for (int k = 0; k < j; ++k) {
                    sum -= cholesky[i][k] * cholesky[j][k];
                }
                if (i == j) {
                    if (!(sum > 0.0) || !Double.isFinite(sum)) {
                        throw new IllegalArgumentException("Matrix is not positive definite");
                    }
                    final double diagonal = Math.sqrt(sum);
                    cholesky[i][j] = diagonal;
                    logDet += Math.log(diagonal);
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
                double sum = (row == column) ? 1.0 : 0.0;
                for (int k = 0; k < row; ++k) {
                    sum -= cholesky[row][k] * lowerInverse[k][column];
                }
                lowerInverse[row][column] = sum / cholesky[row][row];
            }
        }

        transposeSquare(lowerInverse, lowerInverseTranspose, dimension);
        multiplySquare(lowerInverseTranspose, lowerInverse, inverseOut, dimension);
        symmetrize(inverseOut, dimension);
        return logDet;
    }

    private static void validateSquare(final double[][] matrix,
                                       final int dimension,
                                       final String label) {
        if (matrix.length != dimension) {
            throw new IllegalArgumentException(label + " row dimension mismatch");
        }
        for (double[] row : matrix) {
            if (row.length != dimension) {
                throw new IllegalArgumentException(label + " must be square");
            }
        }
    }

    private static void copy(final double[][] source, final double[][] target) {
        for (int i = 0; i < source.length; ++i) {
            System.arraycopy(source[i], 0, target[i], 0, source[i].length);
        }
    }

    private static void copySquare(final double[][] source,
                                   final double[][] target,
                                   final int dimension) {
        for (int i = 0; i < dimension; ++i) {
            System.arraycopy(source[i], 0, target[i], 0, dimension);
        }
    }

    private static void symmetrizeCopy(final double[][] source,
                                       final double[][] target,
                                       final int dimension) {
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                target[i][j] = 0.5 * (source[i][j] + source[j][i]);
            }
        }
    }

    private static void transpose(final double[][] source, final double[][] target) {
        for (int i = 0; i < source.length; ++i) {
            for (int j = 0; j < source[i].length; ++j) {
                target[j][i] = source[i][j];
            }
        }
    }

    private static void transposeSquare(final double[][] source,
                                        final double[][] target,
                                        final int dimension) {
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                target[j][i] = source[i][j];
            }
        }
    }

    private static void multiply(final double[][] left,
                                 final double[][] right,
                                 final double[][] out) {
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

    private static void multiplySquare(final double[][] left,
                                       final double[][] right,
                                       final double[][] out,
                                       final int dimension) {
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                double sum = 0.0;
                for (int k = 0; k < dimension; ++k) {
                    sum += left[i][k] * right[k][j];
                }
                out[i][j] = sum;
            }
        }
    }

    private static void multiply(final double[][] matrix,
                                 final double[] vector,
                                 final double[] out) {
        for (int i = 0; i < matrix.length; ++i) {
            double sum = 0.0;
            for (int j = 0; j < vector.length; ++j) {
                sum += matrix[i][j] * vector[j];
            }
            out[i] = sum;
        }
    }

    private static double dot(final double[] left, final double[] right) {
        double sum = 0.0;
        for (int i = 0; i < left.length; ++i) {
            sum += left[i] * right[i];
        }
        return sum;
    }

    private static void scaleInPlace(final double[][] matrix, final double factor) {
        for (int i = 0; i < matrix.length; ++i) {
            for (int j = 0; j < matrix[i].length; ++j) {
                matrix[i][j] *= factor;
            }
        }
    }

    private static void scaleInPlace(final double[] vector, final double factor) {
        for (int i = 0; i < vector.length; ++i) {
            vector[i] *= factor;
        }
    }

    private static void symmetrize(final double[][] matrix, final int dimension) {
        for (int i = 0; i < dimension; ++i) {
            for (int j = i + 1; j < dimension; ++j) {
                final double average = 0.5 * (matrix[i][j] + matrix[j][i]);
                matrix[i][j] = average;
                matrix[j][i] = average;
            }
        }
    }

    private static void addDiagonalJitter(final double[][] matrix,
                                          final double jitter,
                                          final int dimension) {
        for (int i = 0; i < dimension; ++i) {
            matrix[i][i] += jitter;
        }
    }

    private static double maxAbsDiagonal(final double[][] matrix, final int dimension) {
        double max = 0.0;
        for (int i = 0; i < dimension; ++i) {
            max = Math.max(max, Math.abs(matrix[i][i]));
        }
        return max;
    }

    private static double gershgorinLowerBound(final double[][] matrix, final int dimension) {
        double lower = Double.POSITIVE_INFINITY;
        for (int i = 0; i < dimension; ++i) {
            double radius = 0.0;
            for (int j = 0; j < dimension; ++j) {
                if (j != i) {
                    radius += Math.abs(matrix[i][j]);
                }
            }
            lower = Math.min(lower, matrix[i][i] - radius);
        }
        return lower;
    }

    private static Workspace workspaceFor(final int dimension) {
        final Workspace workspace = WORKSPACE.get();
        workspace.ensureDimension(dimension);
        return workspace;
    }

    private static final class Workspace {
        double[][] cholesky = new double[0][0];
        double[][] lowerInverse = new double[0][0];
        double[][] lowerInverseTranspose = new double[0][0];
        double[][] symmetric = new double[0][0];
        double[][] adjusted = new double[0][0];

        void ensureDimension(final int dimension) {
            if (cholesky.length >= dimension) {
                return;
            }
            cholesky = new double[dimension][dimension];
            lowerInverse = new double[dimension][dimension];
            lowerInverseTranspose = new double[dimension][dimension];
            symmetric = new double[dimension][dimension];
            adjusted = new double[dimension][dimension];
        }
    }
}
