package dr.evomodel.treedatalikelihood.continuous.canonical.message;

/**
 * Utilities for converting expectation-form Gaussian quantities into canonical form.
 *
 * <p>Public API still accepts {@code double[][]} covariance matrices (these come from
 * the moment-form world outside the canonical path). Internal computation uses flat
 * row-major arrays for efficiency.
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

        final Workspace workspace = workspaceFor(dimension);
        final double[] flatCov = workspace.inputFlat;
        copyToFlat(covariance, flatCov, dimension);

        final double logDet = invertSymmetricPositiveDefinite(flatCov, out.precision, dimension, workspace);
        multiplyMV(out.precision, mean, out.information, dimension);
        out.logNormalizer =
                0.5 * (dimension * Math.log(2.0 * Math.PI) + logDet + dot(mean, out.information, dimension));
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

        final Workspace workspace = workspaceFor(dimension);
        final double[] flatCov = workspace.inputFlat;
        final double[] precision = workspace.cholesky; // reuse scratch slot
        copyToFlat(transitionCovariance, flatCov, dimension);

        final double logDet = invertSymmetricPositiveDefinite(flatCov, precision, dimension, workspace);

        final double[] transitionTranspose = workspace.lowerInverse;
        final double[] tmp = workspace.lowerInverseTranspose;

        transposeFromRagged(transitionMatrix, transitionTranspose, dimension);
        multiplyMM(transitionTranspose, precision, tmp, dimension);
        multiplyMMFromRagged(tmp, transitionMatrix, out.precisionXX, dimension);

        multiplyMM(transitionTranspose, precision, out.precisionXY, dimension);
        scaleFlat(out.precisionXY, -1.0, dimension * dimension);

        multiplyMMFromRagged(precision, transitionMatrix, out.precisionYX, dimension);
        scaleFlat(out.precisionYX, -1.0, dimension * dimension);

        System.arraycopy(precision, 0, out.precisionYY, 0, dimension * dimension);

        multiplyMVFromRagged(precision, transitionOffset, out.informationY, dimension);
        multiplyMV(transitionTranspose, out.informationY, out.informationX, dimension);
        scaleVec(out.informationX, -1.0, dimension);

        out.logNormalizer =
                0.5 * (dimension * Math.log(2.0 * Math.PI)
                        + logDet
                        + dot(transitionOffset, out.informationY, dimension));
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
        final double[] flatCov = workspace.inputFlat;
        try {
            invertSymmetricPositiveDefiniteNoLogDet(state.precision, flatCov, dimension, workspace);
        } catch (IllegalArgumentException ignored) {
            invertSymmetricPositiveDefiniteWithJitterNoLogDet(state.precision, flatCov, dimension, workspace);
        }
        copyFromFlat(flatCov, covarianceOut, dimension);
        multiplyMV(flatCov, state.information, meanOut, dimension);
    }

    private static void invertSymmetricPositiveDefiniteWithJitter(final double[] matrix,
                                                                   final double[] inverseOut,
                                                                   final int dimension,
                                                                   final Workspace workspace) {
        final double[] symmetric = workspace.symmetric;
        final double[] adjusted  = workspace.adjusted;
        symmetrizeCopyFlat(matrix, symmetric, dimension);

        final double jitterBase = Math.max(
                SYMMETRIC_JITTER_ABSOLUTE,
                SYMMETRIC_JITTER_RELATIVE * Math.max(1.0, maxAbsDiagonal(symmetric, dimension)));

        double jitter = 0.0;
        double lowerBound = Double.NaN;
        for (int attempt = 0; attempt < 12; ++attempt) {
            System.arraycopy(symmetric, 0, adjusted, 0, dimension * dimension);
            if (jitter > 0.0) addDiagonalJitter(adjusted, jitter, dimension);
            try {
                invertSymmetricPositiveDefinite(adjusted, inverseOut, dimension, workspace);
                return;
            } catch (IllegalArgumentException ignored) {
                // retry
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

    private static void invertSymmetricPositiveDefiniteWithJitterNoLogDet(final double[] matrix,
                                                                          final double[] inverseOut,
                                                                          final int dimension,
                                                                          final Workspace workspace) {
        final double[] symmetric = workspace.symmetric;
        final double[] adjusted  = workspace.adjusted;
        symmetrizeCopyFlat(matrix, symmetric, dimension);

        final double jitterBase = Math.max(
                SYMMETRIC_JITTER_ABSOLUTE,
                SYMMETRIC_JITTER_RELATIVE * Math.max(1.0, maxAbsDiagonal(symmetric, dimension)));

        double jitter = 0.0;
        double lowerBound = Double.NaN;
        for (int attempt = 0; attempt < 12; ++attempt) {
            System.arraycopy(symmetric, 0, adjusted, 0, dimension * dimension);
            if (jitter > 0.0) addDiagonalJitter(adjusted, jitter, dimension);
            try {
                invertSymmetricPositiveDefiniteNoLogDet(adjusted, inverseOut, dimension, workspace);
                return;
            } catch (IllegalArgumentException ignored) {
                // retry
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

    private static double invertSymmetricPositiveDefinite(final double[] matrix,
                                                          final double[] inverseOut,
                                                          final int dimension,
                                                          final Workspace workspace) {
        return invertSymmetricPositiveDefiniteInternal(matrix, inverseOut, dimension, workspace, true);
    }

    private static void invertSymmetricPositiveDefiniteNoLogDet(final double[] matrix,
                                                                final double[] inverseOut,
                                                                final int dimension,
                                                                final Workspace workspace) {
        invertSymmetricPositiveDefiniteInternal(matrix, inverseOut, dimension, workspace, false);
    }

    private static double invertSymmetricPositiveDefiniteInternal(final double[] matrix,
                                                                  final double[] inverseOut,
                                                                  final int dimension,
                                                                  final Workspace workspace,
                                                                  final boolean computeLogDet) {
        final double[] cholesky              = workspace.cholesky;
        final double[] lowerInverse          = workspace.lowerInverse;
        final double[] lowerInverseTranspose = workspace.lowerInverseTranspose;
        System.arraycopy(matrix, 0, cholesky, 0, dimension * dimension);
        double logDet = 0.0;

        for (int i = 0; i < dimension; ++i) {
            final int iOff = i * dimension;
            for (int j = 0; j <= i; ++j) {
                double sum = cholesky[iOff + j];
                for (int k = 0; k < j; ++k) {
                    sum -= cholesky[iOff + k] * cholesky[j * dimension + k];
                }
                if (i == j) {
                    if (!(sum > 0.0) || !Double.isFinite(sum)) {
                        throw new IllegalArgumentException("Matrix is not positive definite");
                    }
                    final double diagonal = Math.sqrt(sum);
                    cholesky[iOff + j] = diagonal;
                    if (computeLogDet) {
                        logDet += Math.log(diagonal);
                    }
                } else {
                    cholesky[iOff + j] = sum / cholesky[j * dimension + j];
                }
            }
            for (int j = i + 1; j < dimension; ++j) {
                cholesky[iOff + j] = 0.0;
            }
        }
        if (computeLogDet) {
            logDet *= 2.0;
        }

        for (int column = 0; column < dimension; ++column) {
            for (int row = 0; row < dimension; ++row) {
                double sum = (row == column) ? 1.0 : 0.0;
                for (int k = 0; k < row; ++k) {
                    sum -= cholesky[row * dimension + k] * lowerInverse[k * dimension + column];
                }
                lowerInverse[row * dimension + column] = sum / cholesky[row * dimension + row];
            }
        }

        transposeSquareFlat(lowerInverse, lowerInverseTranspose, dimension);
        multiplyMM(lowerInverseTranspose, lowerInverse, inverseOut, dimension);
        symmetrizeFlat(inverseOut, dimension);
        return logDet;
    }

    // -----------------------------------------------------------------------
    // Private flat-array helpers
    // -----------------------------------------------------------------------

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

    private static void copyToFlat(final double[][] source, final double[] target, final int dimension) {
        for (int i = 0; i < dimension; ++i) {
            System.arraycopy(source[i], 0, target, i * dimension, dimension);
        }
    }

    private static void copyFromFlat(final double[] source, final double[][] target, final int dimension) {
        for (int i = 0; i < dimension; ++i) {
            System.arraycopy(source, i * dimension, target[i], 0, dimension);
        }
    }

    private static void symmetrizeCopyFlat(final double[] source,
                                           final double[] target,
                                           final int dimension) {
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                target[i * dimension + j] = 0.5 * (source[i * dimension + j] + source[j * dimension + i]);
            }
        }
    }

    private static void transposeFromRagged(final double[][] source,
                                            final double[] target,
                                            final int dimension) {
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                target[j * dimension + i] = source[i][j];
            }
        }
    }

    private static void transposeSquareFlat(final double[] source,
                                            final double[] target,
                                            final int dimension) {
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                target[j * dimension + i] = source[i * dimension + j];
            }
        }
    }

    private static void multiplyMM(final double[] left,
                                   final double[] right,
                                   final double[] out,
                                   final int dimension) {
        for (int i = 0; i < dimension; ++i) {
            final int iOff = i * dimension;
            for (int j = 0; j < dimension; ++j) {
                double sum = 0.0;
                for (int k = 0; k < dimension; ++k) {
                    sum += left[iOff + k] * right[k * dimension + j];
                }
                out[iOff + j] = sum;
            }
        }
    }

    /** Multiply flat × ragged-matrix → flat. */
    private static void multiplyMMFromRagged(final double[] left,
                                             final double[][] right,
                                             final double[] out,
                                             final int dimension) {
        for (int i = 0; i < dimension; ++i) {
            final int iOff = i * dimension;
            for (int j = 0; j < dimension; ++j) {
                double sum = 0.0;
                for (int k = 0; k < dimension; ++k) {
                    sum += left[iOff + k] * right[k][j];
                }
                out[iOff + j] = sum;
            }
        }
    }

    private static void multiplyMV(final double[] matrix,
                                   final double[] vector,
                                   final double[] out,
                                   final int dimension) {
        for (int i = 0; i < dimension; ++i) {
            final int iOff = i * dimension;
            double sum = 0.0;
            for (int j = 0; j < dimension; ++j) {
                sum += matrix[iOff + j] * vector[j];
            }
            out[i] = sum;
        }
    }

    /** Multiply ragged-matrix × vector → flat output. */
    private static void multiplyMVFromRagged(final double[] matrix,
                                              final double[] vector,
                                              final double[] out,
                                              final int dimension) {
        // matrix is already flat here
        multiplyMV(matrix, vector, out, dimension);
    }

    private static double dot(final double[] left, final double[] right, final int dimension) {
        double sum = 0.0;
        for (int i = 0; i < dimension; ++i) {
            sum += left[i] * right[i];
        }
        return sum;
    }

    private static void scaleFlat(final double[] array, final double factor, final int length) {
        for (int k = 0; k < length; ++k) {
            array[k] *= factor;
        }
    }

    private static void scaleVec(final double[] vector, final double factor, final int length) {
        for (int k = 0; k < length; ++k) {
            vector[k] *= factor;
        }
    }

    private static void symmetrizeFlat(final double[] matrix, final int dimension) {
        for (int i = 0; i < dimension; ++i) {
            for (int j = i + 1; j < dimension; ++j) {
                final double avg = 0.5 * (matrix[i * dimension + j] + matrix[j * dimension + i]);
                matrix[i * dimension + j] = avg;
                matrix[j * dimension + i] = avg;
            }
        }
    }

    private static void addDiagonalJitter(final double[] matrix,
                                          final double jitter,
                                          final int dimension) {
        for (int i = 0; i < dimension; ++i) {
            matrix[i * dimension + i] += jitter;
        }
    }

    private static double maxAbsDiagonal(final double[] matrix, final int dimension) {
        double max = 0.0;
        for (int i = 0; i < dimension; ++i) {
            max = Math.max(max, Math.abs(matrix[i * dimension + i]));
        }
        return max;
    }

    private static double gershgorinLowerBound(final double[] matrix, final int dimension) {
        double lower = Double.POSITIVE_INFINITY;
        for (int i = 0; i < dimension; ++i) {
            final int iOff = i * dimension;
            double radius = 0.0;
            for (int j = 0; j < dimension; ++j) {
                if (j != i) radius += Math.abs(matrix[iOff + j]);
            }
            lower = Math.min(lower, matrix[iOff + i] - radius);
        }
        return lower;
    }

    private static Workspace workspaceFor(final int dimension) {
        final Workspace workspace = WORKSPACE.get();
        workspace.ensureDimension(dimension);
        return workspace;
    }

    private static final class Workspace {
        double[] inputFlat          = new double[0];
        double[] cholesky           = new double[0];
        double[] lowerInverse       = new double[0];
        double[] lowerInverseTranspose = new double[0];
        double[] symmetric          = new double[0];
        double[] adjusted           = new double[0];

        void ensureDimension(final int dimension) {
            final int d2 = dimension * dimension;
            if (cholesky.length >= d2) return;
            inputFlat              = new double[d2];
            cholesky               = new double[d2];
            lowerInverse           = new double[d2];
            lowerInverseTranspose  = new double[d2];
            symmetric              = new double[d2];
            adjusted               = new double[d2];
        }
    }
}
