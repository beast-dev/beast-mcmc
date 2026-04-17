package dr.inference.timeseries.engine.gaussian;

import dr.inference.timeseries.representation.CanonicalGaussianState;
import dr.inference.timeseries.representation.CanonicalGaussianTransition;

/**
 * Shared canonical-form Gaussian message operations.
 *
 * <p>This utility centralizes the Schur-complement algebra used by canonical
 * filtering, smoothing, and tree message passing. It operates directly on
 * {@link CanonicalGaussianState} and {@link CanonicalGaussianTransition}.
 */
public final class CanonicalGaussianMessageOps {

    private static final double SYMMETRIC_JITTER_RELATIVE = 1.0e-12;
    private static final double SYMMETRIC_JITTER_ABSOLUTE = 1.0e-12;

    public static final class Workspace {
        private final int dimension;
        final double[][] matrix1;
        final double[][] matrix2;
        final double[][] matrix3;
        final double[][] matrix4;
        final double[] vector1;
        final double[] vector2;

        public Workspace(final int dimension) {
            if (dimension < 1) {
                throw new IllegalArgumentException("dimension must be positive");
            }
            this.dimension = dimension;
            this.matrix1 = new double[dimension][dimension];
            this.matrix2 = new double[dimension][dimension];
            this.matrix3 = new double[dimension][dimension];
            this.matrix4 = new double[dimension][dimension];
            this.vector1 = new double[dimension];
            this.vector2 = new double[dimension];
        }

        public int getDimension() {
            return dimension;
        }
    }

    private CanonicalGaussianMessageOps() { }

    public static void clearState(final CanonicalGaussianState state) {
        final int dim = state.getDimension();
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                state.precision[i][j] = 0.0;
            }
            state.information[i] = 0.0;
        }
        state.logNormalizer = 0.0;
    }

    public static void copyState(final CanonicalGaussianState source, final CanonicalGaussianState target) {
        final int dim = source.getDimension();
        for (int i = 0; i < dim; i++) {
            System.arraycopy(source.precision[i], 0, target.precision[i], 0, dim);
            target.information[i] = source.information[i];
        }
        target.logNormalizer = source.logNormalizer;
    }

    public static void combineStates(final CanonicalGaussianState left,
                                     final CanonicalGaussianState right,
                                     final CanonicalGaussianState out) {
        final int dim = left.getDimension();
        for (int i = 0; i < dim; i++) {
            out.information[i] = left.information[i] + right.information[i];
            for (int j = 0; j < dim; j++) {
                out.precision[i][j] = left.precision[i][j] + right.precision[i][j];
            }
        }
        out.logNormalizer = left.logNormalizer + right.logNormalizer;
    }

    public static void combineStateInPlace(final CanonicalGaussianState accumulator,
                                           final CanonicalGaussianState increment) {
        final int dim = accumulator.getDimension();
        for (int i = 0; i < dim; i++) {
            accumulator.information[i] += increment.information[i];
            for (int j = 0; j < dim; j++) {
                accumulator.precision[i][j] += increment.precision[i][j];
            }
        }
        accumulator.logNormalizer += increment.logNormalizer;
    }

    public static void pushForward(final CanonicalGaussianState previous,
                                   final CanonicalGaussianTransition transition,
                                   final Workspace workspace,
                                   final CanonicalGaussianState out) {
        final int d = previous.getDimension();
        ensureDimension(workspace, d);

        final double[][] a = workspace.matrix1;
        final double[][] aInverse = workspace.matrix2;
        final double[][] temp = workspace.matrix3;
        final double[][] temp2 = workspace.matrix4;
        final double[] h = workspace.vector1;
        final double[] tempv = workspace.vector2;

        for (int i = 0; i < d; ++i) {
            h[i] = previous.information[i] + transition.informationX[i];
            for (int j = 0; j < d; ++j) {
                a[i][j] = previous.precision[i][j] + transition.precisionXX[i][j];
            }
        }

        final double eliminated = normalizedLogNormalizer(a, h, d, workspace, aInverse, tempv);

        KalmanLikelihoodEngine.multiplyMatrixMatrix(aInverse, transition.precisionXY, temp);
        KalmanLikelihoodEngine.multiplyMatrixMatrix(transition.precisionYX, temp, temp2);

        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                out.precision[i][j] = transition.precisionYY[i][j] - temp2[i][j];
            }
        }
        KalmanLikelihoodEngine.symmetrize(out.precision);

        KalmanLikelihoodEngine.multiplyMatrixVector(aInverse, h, tempv, d, d);
        KalmanLikelihoodEngine.multiplyMatrixVector(transition.precisionYX, tempv, h, d, d);
        for (int i = 0; i < d; ++i) {
            out.information[i] = transition.informationY[i] - h[i];
        }
        out.logNormalizer = previous.logNormalizer + transition.logNormalizer - eliminated;
    }

    public static void pushBackward(final CanonicalGaussianState next,
                                    final CanonicalGaussianTransition transition,
                                    final Workspace workspace,
                                    final CanonicalGaussianState out) {
        final int d = next.getDimension();
        ensureDimension(workspace, d);

        final double[][] a = workspace.matrix1;
        final double[][] aInverse = workspace.matrix2;
        final double[][] temp = workspace.matrix3;
        final double[][] temp2 = workspace.matrix4;
        final double[] h = workspace.vector1;
        final double[] tempv = workspace.vector2;

        for (int i = 0; i < d; ++i) {
            h[i] = next.information[i] + transition.informationY[i];
            for (int j = 0; j < d; ++j) {
                a[i][j] = next.precision[i][j] + transition.precisionYY[i][j];
            }
        }

        final double eliminated = normalizedLogNormalizer(a, h, d, workspace, aInverse, tempv);

        KalmanLikelihoodEngine.multiplyMatrixMatrix(transition.precisionXY, aInverse, temp);
        KalmanLikelihoodEngine.multiplyMatrixMatrix(temp, transition.precisionYX, temp2);

        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                out.precision[i][j] = transition.precisionXX[i][j] - temp2[i][j];
            }
        }
        KalmanLikelihoodEngine.symmetrize(out.precision);

        KalmanLikelihoodEngine.multiplyMatrixVector(aInverse, h, tempv, d, d);
        KalmanLikelihoodEngine.multiplyMatrixVector(transition.precisionXY, tempv, h, d, d);
        for (int i = 0; i < d; ++i) {
            out.information[i] = transition.informationX[i] - h[i];
        }
        out.logNormalizer = next.logNormalizer + transition.logNormalizer - eliminated;
    }

    public static void conditionOnObservedSecondBlock(final CanonicalGaussianTransition transition,
                                                      final double[] observed,
                                                      final CanonicalGaussianState out) {
        final int d = transition.getDimension();
        for (int i = 0; i < d; ++i) {
            double info = transition.informationX[i];
            for (int j = 0; j < d; ++j) {
                out.precision[i][j] = transition.precisionXX[i][j];
                info -= transition.precisionXY[i][j] * observed[j];
            }
            out.information[i] = info;
        }
        out.logNormalizer = transition.logNormalizer
                + 0.5 * quadraticForm(transition.precisionYY, observed, d)
                - dot(transition.informationY, observed, d);
    }

    public static void conditionOnObservedFirstBlock(final CanonicalGaussianTransition transition,
                                                     final double[] observed,
                                                     final CanonicalGaussianState out) {
        final int d = transition.getDimension();
        for (int i = 0; i < d; ++i) {
            double info = transition.informationY[i];
            for (int j = 0; j < d; ++j) {
                out.precision[i][j] = transition.precisionYY[i][j];
                info -= transition.precisionYX[i][j] * observed[j];
            }
            out.information[i] = info;
        }
        out.logNormalizer = transition.logNormalizer
                + 0.5 * quadraticForm(transition.precisionXX, observed, d)
                - dot(transition.informationX, observed, d);
    }

    /**
     * Conditions a canonical transition on an exact observation of a subset of the child
     * coordinates and marginalizes the remaining child coordinates.
     *
     * <p>The output is a canonical Gaussian state over the parent block only. The
     * {@code observedValues} array is indexed in the original child coordinate system; only
     * entries referenced by {@code observedIndices} are read.
     */
    public static void conditionOnPartiallyObservedSecondBlock(final CanonicalGaussianTransition transition,
                                                               final double[] observedValues,
                                                               final int[] observedIndices,
                                                               final int observedCount,
                                                               final int[] missingIndices,
                                                               final int missingCount,
                                                               final Workspace workspace,
                                                               final CanonicalGaussianState out) {
        final int d = transition.getDimension();
        ensureDimension(workspace, d);
        if (observedCount < 0 || missingCount < 0 || observedCount + missingCount != d) {
            throw new IllegalArgumentException("Observed/missing partition must cover the child dimension.");
        }
        if (observedCount == 0) {
            clearState(out);
            return;
        }
        if (missingCount == 0) {
            conditionOnObservedSecondBlock(transition, observedValues, out);
            return;
        }

        final double[][] missingPrecision = workspace.matrix1;
        final double[][] missingPrecisionInverse = workspace.matrix2;
        final double[][] xyTimesMissingInverse = workspace.matrix3;
        final double[] missingInformation = workspace.vector1;
        final double[] missingMean = workspace.vector2;

        for (int i = 0; i < d; ++i) {
            double information = transition.informationX[i];
            for (int j = 0; j < d; ++j) {
                out.precision[i][j] = transition.precisionXX[i][j];
            }
            for (int observed = 0; observed < observedCount; ++observed) {
                final int observedIndex = observedIndices[observed];
                information -= transition.precisionXY[i][observedIndex] * observedValues[observedIndex];
            }
            out.information[i] = information;
        }

        for (int row = 0; row < missingCount; ++row) {
            final int missingRow = missingIndices[row];
            double information = transition.informationY[missingRow];
            for (int col = 0; col < missingCount; ++col) {
                missingPrecision[row][col] = transition.precisionYY[missingRow][missingIndices[col]];
            }
            for (int observed = 0; observed < observedCount; ++observed) {
                final int observedIndex = observedIndices[observed];
                information -= transition.precisionYY[missingRow][observedIndex] * observedValues[observedIndex];
            }
            missingInformation[row] = information;
        }

        final double eliminated;
        try {
            eliminated = normalizedLogNormalizer(
                    missingPrecision,
                    missingInformation,
                    missingCount,
                    workspace,
                    missingPrecisionInverse,
                    missingMean);
        } catch (IllegalArgumentException e) {
            throw new IllegalArgumentException(
                    "Failed canonical partial-tip elimination with missing block "
                            + formatSquare(missingPrecision, missingCount)
                            + " and missing information "
                            + formatVector(missingInformation, missingCount),
                    e);
        }

        for (int i = 0; i < d; ++i) {
            double informationCorrection = 0.0;
            for (int missing = 0; missing < missingCount; ++missing) {
                final int missingIndex = missingIndices[missing];
                informationCorrection += transition.precisionXY[i][missingIndex] * missingMean[missing];
            }
            out.information[i] -= informationCorrection;

            for (int missing = 0; missing < missingCount; ++missing) {
                double sum = 0.0;
                for (int otherMissing = 0; otherMissing < missingCount; ++otherMissing) {
                    sum += transition.precisionXY[i][missingIndices[otherMissing]]
                            * missingPrecisionInverse[otherMissing][missing];
                }
                xyTimesMissingInverse[i][missing] = sum;
            }
        }

        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                double precisionCorrection = 0.0;
                for (int missing = 0; missing < missingCount; ++missing) {
                    precisionCorrection += xyTimesMissingInverse[i][missing]
                            * transition.precisionYX[missingIndices[missing]][j];
                }
                out.precision[i][j] -= precisionCorrection;
            }
        }
        KalmanLikelihoodEngine.symmetrize(out.precision);
        out.logNormalizer = transition.logNormalizer
                + observedQuadraticConstant(transition, observedValues, observedIndices, observedCount)
                - eliminated;
    }

    public static void buildPairPosterior(final CanonicalGaussianState firstState,
                                          final CanonicalGaussianTransition transition,
                                          final CanonicalGaussianState secondState,
                                          final CanonicalGaussianState pairOut) {
        final int d = firstState.getDimension();
        for (int i = 0; i < d; ++i) {
            pairOut.information[i] = firstState.information[i] + transition.informationX[i];
            pairOut.information[d + i] = secondState.information[i] + transition.informationY[i];
            for (int j = 0; j < d; ++j) {
                pairOut.precision[i][j] = firstState.precision[i][j] + transition.precisionXX[i][j];
                pairOut.precision[i][d + j] = transition.precisionXY[i][j];
                pairOut.precision[d + i][j] = transition.precisionYX[i][j];
                pairOut.precision[d + i][d + j] = secondState.precision[i][j] + transition.precisionYY[i][j];
            }
        }
        pairOut.logNormalizer = firstState.logNormalizer + transition.logNormalizer + secondState.logNormalizer;
    }

    public static void marginalizeFirstBlock(final CanonicalGaussianState pairState,
                                             final int blockDimension,
                                             final Workspace workspace,
                                             final CanonicalGaussianState out) {
        marginalize(pairState, blockDimension, true, workspace, out);
    }

    public static void marginalizeSecondBlock(final CanonicalGaussianState pairState,
                                              final int blockDimension,
                                              final Workspace workspace,
                                              final CanonicalGaussianState out) {
        marginalize(pairState, blockDimension, false, workspace, out);
    }

    public static double normalizedLogNormalizer(final CanonicalGaussianState state,
                                                 final Workspace workspace) {
        final int d = state.getDimension();
        ensureDimension(workspace, d);
        return normalizedLogNormalizer(state.precision, state.information, d, workspace, workspace.matrix1, workspace.vector1);
    }

    public static double normalizationShift(final CanonicalGaussianState state,
                                            final Workspace workspace) {
        return normalizedLogNormalizer(state, workspace) - state.logNormalizer;
    }

    private static void marginalize(final CanonicalGaussianState pairState,
                                    final int d,
                                    final boolean keepFirstBlock,
                                    final Workspace workspace,
                                    final CanonicalGaussianState out) {
        ensureDimension(workspace, d);
        final double[][] elimPrecision = workspace.matrix1;
        final double[][] elimPrecisionInverse = workspace.matrix2;
        final double[][] cross1 = workspace.matrix3;
        final double[][] cross2 = workspace.matrix4;
        final double[] elimInformation = workspace.vector1;
        final double[] tempv = workspace.vector2;

        if (keepFirstBlock) {
            fillLowerRightBlock(pairState.precision, d, elimPrecision);
            fillUpperRightBlock(pairState.precision, d, cross1);
            fillLowerLeftBlock(pairState.precision, d, cross2);
            fillSecondBlock(pairState.information, d, elimInformation);
        } else {
            fillUpperLeftBlock(pairState.precision, d, elimPrecision);
            fillLowerLeftBlock(pairState.precision, d, cross1);
            fillUpperRightBlock(pairState.precision, d, cross2);
            fillFirstBlock(pairState.information, d, elimInformation);
        }

        final double eliminated = normalizedLogNormalizer(elimPrecision, elimInformation, d,
                workspace, elimPrecisionInverse, tempv);

        KalmanLikelihoodEngine.multiplyMatrixMatrix(cross1, elimPrecisionInverse, workspace.matrix3);
        KalmanLikelihoodEngine.multiplyMatrixMatrix(workspace.matrix3, cross2, workspace.matrix4);

        if (keepFirstBlock) {
            fillUpperLeftBlock(pairState.precision, d, out.precision);
            fillFirstBlock(pairState.information, d, out.information);
        } else {
            fillLowerRightBlock(pairState.precision, d, out.precision);
            fillSecondBlock(pairState.information, d, out.information);
        }

        subtractMatrixInPlace(out.precision, workspace.matrix4, d);
        KalmanLikelihoodEngine.symmetrize(out.precision);

        KalmanLikelihoodEngine.multiplyMatrixVector(elimPrecisionInverse, elimInformation, tempv, d, d);
        KalmanLikelihoodEngine.multiplyMatrixVector(cross1, tempv, elimInformation, d, d);
        for (int i = 0; i < d; ++i) {
            out.information[i] -= elimInformation[i];
        }

        out.logNormalizer = pairState.logNormalizer - eliminated;
    }

    private static double normalizedLogNormalizer(final double[][] precision,
                                                  final double[] information,
                                                  final int dimension,
                                                  final Workspace workspace,
                                                  final double[][] inverseOut,
                                                  final double[] tempVector) {
        final double logDet = invertPositiveDefinite(precision, inverseOut, dimension, workspace);
        KalmanLikelihoodEngine.multiplyMatrixVector(inverseOut, information, tempVector, dimension, dimension);
        final double quadratic = dot(information, tempVector, dimension);
        return 0.5 * (dimension * KalmanLikelihoodEngine.LOG_TWO_PI - logDet + quadratic);
    }

    private static double invertPositiveDefinite(final double[][] matrix,
                                                 final double[][] inverseOut,
                                                 final int dimension,
                                                 final Workspace workspace) {
        KalmanLikelihoodEngine.copyMatrix(matrix, workspace.matrix4, dimension, dimension);
        symmetrizeSquare(workspace.matrix4, dimension);
        try {
            return invertPositiveDefiniteFromSymmetricCopy(workspace.matrix4, inverseOut, dimension);
        } catch (IllegalArgumentException ignored) {
            // Retry with progressively larger diagonal jitter when numerical symmetry/PD issues arise.
        }

        KalmanLikelihoodEngine.copyMatrix(workspace.matrix4, workspace.matrix3, dimension, dimension);
        final double jitterBase = Math.max(
                SYMMETRIC_JITTER_ABSOLUTE,
                SYMMETRIC_JITTER_RELATIVE * Math.max(1.0, maxAbsDiagonal(workspace.matrix3, dimension)));

        double jitter = 0.0;
        double lowerBound = Double.NaN;
        for (int attempt = 0; attempt < 12; ++attempt) {
            KalmanLikelihoodEngine.copyMatrix(workspace.matrix3, workspace.matrix4, dimension, dimension);
            if (jitter > 0.0) {
                addDiagonalJitter(workspace.matrix4, jitter, dimension);
            }
            try {
                return invertPositiveDefiniteFromSymmetricCopy(workspace.matrix4, inverseOut, dimension);
            } catch (IllegalArgumentException ignored) {
                // Increase jitter and retry.
            }

            if (jitter == 0.0) {
                if (Double.isNaN(lowerBound)) {
                    lowerBound = gershgorinLowerBound(workspace.matrix3, dimension);
                }
                jitter = lowerBound > 0.0 ? jitterBase : (-lowerBound + jitterBase);
            } else {
                jitter *= 10.0;
            }
        }

        throw new IllegalArgumentException("Matrix is not positive definite (failed robust inversion retries)");
    }

    private static double invertPositiveDefiniteFromSymmetricCopy(final double[][] matrix,
                                                                  final double[][] inverseOut,
                                                                  final int dimension) {
        final KalmanLikelihoodEngine.CholeskyFactor chol = KalmanLikelihoodEngine.cholesky(matrix);
        KalmanLikelihoodEngine.copyMatrix(matrix, inverseOut, dimension, dimension);
        KalmanLikelihoodEngine.invertPositiveDefiniteFromCholesky(inverseOut, chol);
        return chol.logDeterminant();
    }

    private static double observedQuadraticConstant(final CanonicalGaussianTransition transition,
                                                    final double[] observedValues,
                                                    final int[] observedIndices,
                                                    final int observedCount) {
        double quadratic = 0.0;
        double linear = 0.0;
        for (int i = 0; i < observedCount; ++i) {
            final int observedIndexI = observedIndices[i];
            final double observedValueI = observedValues[observedIndexI];
            linear += transition.informationY[observedIndexI] * observedValueI;
            for (int j = 0; j < observedCount; ++j) {
                final int observedIndexJ = observedIndices[j];
                quadratic += observedValueI
                        * transition.precisionYY[observedIndexI][observedIndexJ]
                        * observedValues[observedIndexJ];
            }
        }
        return 0.5 * quadratic - linear;
    }

    private static double quadraticForm(final double[][] matrix, final double[] vector, final int dim) {
        double result = 0.0;
        for (int i = 0; i < dim; i++) {
            double row = 0.0;
            for (int j = 0; j < dim; j++) {
                row += matrix[i][j] * vector[j];
            }
            result += vector[i] * row;
        }
        return result;
    }

    private static double dot(final double[] left, final double[] right, final int dim) {
        double sum = 0.0;
        for (int i = 0; i < dim; i++) {
            sum += left[i] * right[i];
        }
        return sum;
    }

    private static void subtractMatrixInPlace(final double[][] target, final double[][] delta, final int dim) {
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                target[i][j] -= delta[i][j];
            }
        }
    }

    private static void fillUpperLeftBlock(final double[][] source, final int d, final double[][] out) {
        for (int i = 0; i < d; ++i) {
            System.arraycopy(source[i], 0, out[i], 0, d);
        }
    }

    private static void fillUpperRightBlock(final double[][] source, final int d, final double[][] out) {
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                out[i][j] = source[i][d + j];
            }
        }
    }

    private static void fillLowerLeftBlock(final double[][] source, final int d, final double[][] out) {
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                out[i][j] = source[d + i][j];
            }
        }
    }

    private static void fillLowerRightBlock(final double[][] source, final int d, final double[][] out) {
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                out[i][j] = source[d + i][d + j];
            }
        }
    }

    private static void fillFirstBlock(final double[] source, final int d, final double[] out) {
        System.arraycopy(source, 0, out, 0, d);
    }

    private static void fillSecondBlock(final double[] source, final int d, final double[] out) {
        System.arraycopy(source, d, out, 0, d);
    }

    private static void ensureDimension(final Workspace workspace, final int dimension) {
        if (workspace.getDimension() != dimension) {
            throw new IllegalArgumentException(
                    "Workspace dimension mismatch: " + workspace.getDimension() + " vs " + dimension);
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

    private static void symmetrizeSquare(final double[][] matrix, final int dimension) {
        for (int i = 0; i < dimension; ++i) {
            for (int j = i + 1; j < dimension; ++j) {
                final double average = 0.5 * (matrix[i][j] + matrix[j][i]);
                matrix[i][j] = average;
                matrix[j][i] = average;
            }
        }
    }

    private static String formatSquare(final double[][] matrix, final int dimension) {
        final StringBuilder builder = new StringBuilder("[");
        for (int i = 0; i < dimension; ++i) {
            if (i > 0) {
                builder.append(", ");
            }
            builder.append('[');
            for (int j = 0; j < dimension; ++j) {
                if (j > 0) {
                    builder.append(", ");
                }
                builder.append(matrix[i][j]);
            }
            builder.append(']');
        }
        builder.append(']');
        return builder.toString();
    }

    private static String formatVector(final double[] vector, final int dimension) {
        final StringBuilder builder = new StringBuilder("[");
        for (int i = 0; i < dimension; ++i) {
            if (i > 0) {
                builder.append(", ");
            }
            builder.append(vector[i]);
        }
        builder.append(']');
        return builder.toString();
    }
}
