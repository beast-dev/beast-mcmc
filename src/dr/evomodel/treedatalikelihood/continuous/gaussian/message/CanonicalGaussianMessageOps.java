package dr.evomodel.treedatalikelihood.continuous.gaussian.message;

import dr.evomodel.treedatalikelihood.continuous.gaussian.CanonicalGaussianState;
import dr.evomodel.treedatalikelihood.continuous.gaussian.CanonicalGaussianTransition;

/**
 * Shared canonical-form Gaussian message operations.
 *
 * <p>This utility centralizes the Schur-complement algebra used by canonical
 * filtering, smoothing, and tree message passing. It operates directly on
 * {@link CanonicalGaussianState} and {@link CanonicalGaussianTransition}.
 *
 * <p>All matrix fields use flat row-major storage: element {@code (i, j)} of
 * a {@code dim × dim} matrix is at index {@code i * dim + j}.
 */
public final class CanonicalGaussianMessageOps {

    private static final double SYMMETRIC_JITTER_RELATIVE = 1.0e-12;
    private static final double SYMMETRIC_JITTER_ABSOLUTE = 1.0e-12;

    public static final class Workspace {
        private final int dimension;
        final double[] matrix1;
        final double[] matrix2;
        final double[] matrix3;
        final double[] matrix4;
        final double[] matrix5;
        final double[] vector1;
        final double[] vector2;
        final double[] vector3;

        public Workspace(final int dimension) {
            if (dimension < 1) {
                throw new IllegalArgumentException("dimension must be positive");
            }
            this.dimension = dimension;
            final int d2 = dimension * dimension;
            this.matrix1 = new double[d2];
            this.matrix2 = new double[d2];
            this.matrix3 = new double[d2];
            this.matrix4 = new double[d2];
            this.matrix5 = new double[d2];
            this.vector1 = new double[dimension];
            this.vector2 = new double[dimension];
            this.vector3 = new double[dimension];
        }

        public int getDimension() {
            return dimension;
        }
    }

    private CanonicalGaussianMessageOps() { }

    public static void clearState(final CanonicalGaussianState state) {
        final int dim = state.getDimension();
        java.util.Arrays.fill(state.precision, 0, dim * dim, 0.0);
        java.util.Arrays.fill(state.information, 0.0);
        state.logNormalizer = 0.0;
    }

    public static void copyState(final CanonicalGaussianState source, final CanonicalGaussianState target) {
        final int dim = source.getDimension();
        System.arraycopy(source.precision, 0, target.precision, 0, dim * dim);
        System.arraycopy(source.information, 0, target.information, 0, dim);
        target.logNormalizer = source.logNormalizer;
    }

    public static void combineStates(final CanonicalGaussianState left,
                                     final CanonicalGaussianState right,
                                     final CanonicalGaussianState out) {
        final int dim = left.getDimension();
        final int d2 = dim * dim;
        for (int i = 0; i < dim; i++) {
            out.information[i] = left.information[i] + right.information[i];
        }
        for (int k = 0; k < d2; k++) {
            out.precision[k] = left.precision[k] + right.precision[k];
        }
        out.logNormalizer = left.logNormalizer + right.logNormalizer;
    }

    public static void combineStateInPlace(final CanonicalGaussianState accumulator,
                                           final CanonicalGaussianState increment) {
        final int dim = accumulator.getDimension();
        final int d2 = dim * dim;
        for (int i = 0; i < dim; i++) {
            accumulator.information[i] += increment.information[i];
        }
        for (int k = 0; k < d2; k++) {
            accumulator.precision[k] += increment.precision[k];
        }
        accumulator.logNormalizer += increment.logNormalizer;
    }

    public static void pushForward(final CanonicalGaussianState previous,
                                   final CanonicalGaussianTransition transition,
                                   final Workspace workspace,
                                   final CanonicalGaussianState out) {
        final int d = previous.getDimension();
        ensureDimension(workspace, d);

        final double[] a      = workspace.matrix1;
        final double[] aInv   = workspace.matrix2;
        final double[] temp   = workspace.matrix3;
        final double[] temp2  = workspace.matrix4;
        final double[] h      = workspace.vector1;
        final double[] tempv  = workspace.vector2;

        for (int i = 0; i < d; ++i) {
            h[i] = previous.information[i] + transition.informationX[i];
            final int iOff = i * d;
            for (int j = 0; j < d; ++j) {
                a[iOff + j] = previous.precision[iOff + j] + transition.precisionXX[iOff + j];
            }
        }

        final double eliminated = normalizedLogNormalizer(a, h, d, workspace, aInv, tempv);

        GaussianMatrixOps.multiplyMatrixMatrixFlat(aInv, transition.precisionXY, temp, d);
        GaussianMatrixOps.multiplyMatrixMatrixFlat(transition.precisionYX, temp, temp2, d);

        final int d2 = d * d;
        for (int k = 0; k < d2; ++k) {
            out.precision[k] = transition.precisionYY[k] - temp2[k];
        }
        GaussianMatrixOps.symmetrizeFlat(out.precision, d);

        GaussianMatrixOps.multiplyMatrixVectorFlat(aInv, h, tempv, d);
        GaussianMatrixOps.multiplyMatrixVectorFlat(transition.precisionYX, tempv, h, d);
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

        final double[] a      = workspace.matrix1;
        final double[] aInv   = workspace.matrix2;
        final double[] temp   = workspace.matrix3;
        final double[] temp2  = workspace.matrix4;
        final double[] h      = workspace.vector1;
        final double[] tempv  = workspace.vector2;

        for (int i = 0; i < d; ++i) {
            h[i] = next.information[i] + transition.informationY[i];
            final int iOff = i * d;
            for (int j = 0; j < d; ++j) {
                a[iOff + j] = next.precision[iOff + j] + transition.precisionYY[iOff + j];
            }
        }

        final double eliminated = normalizedLogNormalizer(a, h, d, workspace, aInv, tempv);

        GaussianMatrixOps.multiplyMatrixMatrixFlat(transition.precisionXY, aInv, temp, d);
        GaussianMatrixOps.multiplyMatrixMatrixFlat(temp, transition.precisionYX, temp2, d);

        final int d2 = d * d;
        for (int k = 0; k < d2; ++k) {
            out.precision[k] = transition.precisionXX[k] - temp2[k];
        }
        GaussianMatrixOps.symmetrizeFlat(out.precision, d);

        GaussianMatrixOps.multiplyMatrixVectorFlat(aInv, h, tempv, d);
        GaussianMatrixOps.multiplyMatrixVectorFlat(transition.precisionXY, tempv, h, d);
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
            final int iOff = i * d;
            double info = transition.informationX[i];
            for (int j = 0; j < d; ++j) {
                out.precision[iOff + j] = transition.precisionXX[iOff + j];
                info -= transition.precisionXY[iOff + j] * observed[j];
            }
            out.information[i] = info;
        }
        out.logNormalizer = transition.logNormalizer
                + 0.5 * quadraticFormFlat(transition.precisionYY, observed, d)
                - dotFlat(transition.informationY, observed, d);
    }

    public static void conditionOnObservedFirstBlock(final CanonicalGaussianTransition transition,
                                                     final double[] observed,
                                                     final CanonicalGaussianState out) {
        final int d = transition.getDimension();
        for (int i = 0; i < d; ++i) {
            final int iOff = i * d;
            double info = transition.informationY[i];
            for (int j = 0; j < d; ++j) {
                out.precision[iOff + j] = transition.precisionYY[iOff + j];
                info -= transition.precisionYX[iOff + j] * observed[j];
            }
            out.information[i] = info;
        }
        out.logNormalizer = transition.logNormalizer
                + 0.5 * quadraticFormFlat(transition.precisionXX, observed, d)
                - dotFlat(transition.informationX, observed, d);
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

        final double[] missingPrecision        = workspace.matrix1;
        final double[] missingPrecisionInverse  = workspace.matrix2;
        final double[] xyTimesMissingInverse    = workspace.matrix3;
        final double[] missingInformation       = workspace.vector1;
        final double[] missingMean              = workspace.vector2;

        for (int i = 0; i < d; ++i) {
            final int iOff = i * d;
            double information = transition.informationX[i];
            for (int j = 0; j < d; ++j) {
                out.precision[iOff + j] = transition.precisionXX[iOff + j];
            }
            for (int obs = 0; obs < observedCount; ++obs) {
                final int obsIdx = observedIndices[obs];
                information -= transition.precisionXY[iOff + obsIdx] * observedValues[obsIdx];
            }
            out.information[i] = information;
        }

        for (int row = 0; row < missingCount; ++row) {
            final int missingRow = missingIndices[row];
            final int missingRowOff = missingRow * d;
            double information = transition.informationY[missingRow];
            for (int col = 0; col < missingCount; ++col) {
                missingPrecision[row * missingCount + col] =
                        transition.precisionYY[missingRowOff + missingIndices[col]];
            }
            for (int obs = 0; obs < observedCount; ++obs) {
                final int obsIdx = observedIndices[obs];
                information -= transition.precisionYY[missingRowOff + obsIdx] * observedValues[obsIdx];
            }
            missingInformation[row] = information;
        }

        final double eliminated;
        try {
            eliminated = normalizedLogNormalizerSubmatrix(
                    missingPrecision,
                    missingInformation,
                    missingCount,
                    workspace,
                    missingPrecisionInverse,
                    missingMean);
        } catch (IllegalArgumentException e) {
            throw new IllegalArgumentException(
                    "Failed canonical partial-tip elimination with missing block "
                            + formatVectorAsSquare(missingPrecision, missingCount)
                            + " and missing information "
                            + formatVector(missingInformation, missingCount),
                    e);
        }

        for (int i = 0; i < d; ++i) {
            final int iOff = i * d;
            double informationCorrection = 0.0;
            for (int missing = 0; missing < missingCount; ++missing) {
                final int missingIndex = missingIndices[missing];
                informationCorrection += transition.precisionXY[iOff + missingIndex] * missingMean[missing];
            }
            out.information[i] -= informationCorrection;

            for (int missing = 0; missing < missingCount; ++missing) {
                double sum = 0.0;
                for (int otherMissing = 0; otherMissing < missingCount; ++otherMissing) {
                    sum += transition.precisionXY[iOff + missingIndices[otherMissing]]
                            * missingPrecisionInverse[otherMissing * missingCount + missing];
                }
                xyTimesMissingInverse[i * missingCount + missing] = sum;
            }
        }

        for (int i = 0; i < d; ++i) {
            final int iOff = i * d;
            for (int j = 0; j < d; ++j) {
                double precisionCorrection = 0.0;
                for (int missing = 0; missing < missingCount; ++missing) {
                    precisionCorrection += xyTimesMissingInverse[i * missingCount + missing]
                            * transition.precisionYX[missingIndices[missing] * d + j];
                }
                out.precision[iOff + j] -= precisionCorrection;
            }
        }
        GaussianMatrixOps.symmetrizeFlat(out.precision, d);
        out.logNormalizer = transition.logNormalizer
                + observedQuadraticConstant(transition, observedValues, observedIndices, observedCount)
                - eliminated;
    }

    public static void buildPairPosterior(final CanonicalGaussianState firstState,
                                          final CanonicalGaussianTransition transition,
                                          final CanonicalGaussianState secondState,
                                          final CanonicalGaussianState pairOut) {
        final int d = firstState.getDimension();
        final int d2 = 2 * d;
        for (int i = 0; i < d; ++i) {
            pairOut.information[i]     = firstState.information[i]  + transition.informationX[i];
            pairOut.information[d + i] = secondState.information[i] + transition.informationY[i];
            final int iOff = i * d;
            for (int j = 0; j < d; ++j) {
                pairOut.precision[i * d2 + j]           = firstState.precision[iOff + j]  + transition.precisionXX[iOff + j];
                pairOut.precision[i * d2 + (d + j)]     = transition.precisionXY[iOff + j];
                pairOut.precision[(d + i) * d2 + j]     = transition.precisionYX[iOff + j];
                pairOut.precision[(d + i) * d2 + (d + j)] = secondState.precision[iOff + j] + transition.precisionYY[iOff + j];
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
        return normalizedLogNormalizer(state.precision, state.information, d, workspace,
                workspace.matrix1, workspace.vector1);
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
        final double[] elimPrecision        = workspace.matrix1;
        final double[] elimPrecisionInverse = workspace.matrix2;
        final double[] cross1               = workspace.matrix3;
        final double[] cross2               = workspace.matrix4;
        final double[] elimInformation      = workspace.vector1;
        final double[] tempv                = workspace.vector2;

        if (keepFirstBlock) {
            fillLowerRightBlockFlat(pairState.precision, d, elimPrecision);
            fillSecondBlockFlat(pairState.information, d, elimInformation);
        } else {
            fillUpperLeftBlockFlat(pairState.precision, d, elimPrecision);
            fillFirstBlockFlat(pairState.information, d, elimInformation);
        }

        final double eliminated = normalizedLogNormalizer(elimPrecision, elimInformation, d,
                workspace, elimPrecisionInverse, tempv);

        if (keepFirstBlock) {
            fillUpperRightBlockFlat(pairState.precision, d, cross1);
            fillLowerLeftBlockFlat(pairState.precision, d, cross2);
        } else {
            fillLowerLeftBlockFlat(pairState.precision, d, cross1);
            fillUpperRightBlockFlat(pairState.precision, d, cross2);
        }

        GaussianMatrixOps.multiplyMatrixMatrixFlat(cross1, elimPrecisionInverse, workspace.matrix5, d);
        GaussianMatrixOps.multiplyMatrixMatrixFlat(workspace.matrix5, cross2, workspace.matrix1, d);

        if (keepFirstBlock) {
            fillUpperLeftBlockFlat(pairState.precision, d, out.precision);
            fillFirstBlockFlat(pairState.information, d, out.information);
        } else {
            fillLowerRightBlockFlat(pairState.precision, d, out.precision);
            fillSecondBlockFlat(pairState.information, d, out.information);
        }

        subtractMatrixInPlaceFlat(out.precision, workspace.matrix1, d);
        GaussianMatrixOps.symmetrizeFlat(out.precision, d);

        GaussianMatrixOps.multiplyMatrixVectorFlat(elimPrecisionInverse, elimInformation, tempv, d);
        GaussianMatrixOps.multiplyMatrixVectorFlat(cross1, tempv, elimInformation, d);
        for (int i = 0; i < d; ++i) {
            out.information[i] -= elimInformation[i];
        }

        out.logNormalizer = pairState.logNormalizer - eliminated;
    }

    private static double normalizedLogNormalizer(final double[] precision,
                                                  final double[] information,
                                                  final int dimension,
                                                  final Workspace workspace,
                                                  final double[] inverseOut,
                                                  final double[] tempVector) {
        final double logDet = invertPositiveDefiniteFlat(precision, inverseOut, dimension, workspace);
        GaussianMatrixOps.multiplyMatrixVectorFlat(inverseOut, information, tempVector, dimension);
        final double quadratic = dotFlat(information, tempVector, dimension);
        return 0.5 * (dimension * GaussianMatrixOps.LOG_TWO_PI - logDet + quadratic);
    }

    /** Same as {@link #normalizedLogNormalizer} but operates on a compact {@code dim x dim} block
     *  stored in the leading portion of the workspace arrays. */
    private static double normalizedLogNormalizerSubmatrix(final double[] precision,
                                                           final double[] information,
                                                           final int dim,
                                                           final Workspace workspace,
                                                           final double[] inverseOut,
                                                           final double[] tempVector) {
        final double logDet = invertPositiveDefiniteFlatSubmatrix(precision, inverseOut, dim, workspace);
        double quadratic = 0.0;
        for (int i = 0; i < dim; ++i) {
            double row = 0.0;
            for (int j = 0; j < dim; ++j) {
                row += inverseOut[i * dim + j] * information[j];
            }
            tempVector[i] = row;
            quadratic += information[i] * row;
        }
        return 0.5 * (dim * GaussianMatrixOps.LOG_TWO_PI - logDet + quadratic);
    }

    private static double invertPositiveDefiniteFlat(final double[] matrix,
                                                      final double[] inverseOut,
                                                      final int dimension,
                                                      final Workspace workspace) {
        GaussianMatrixOps.copyMatrixFlat(matrix, workspace.matrix4, dimension);
        symmetrizeSquareFlat(workspace.matrix4, dimension);
        double logDet = invertPositiveDefiniteFromSymmetricCopyFlat(
                workspace.matrix4, inverseOut, dimension, workspace.matrix5, workspace.matrix4);
        if (!Double.isNaN(logDet)) {
            return logDet;
        }

        GaussianMatrixOps.copyMatrixFlat(workspace.matrix4, workspace.matrix3, dimension);
        final double jitterBase = Math.max(
                SYMMETRIC_JITTER_ABSOLUTE,
                SYMMETRIC_JITTER_RELATIVE * Math.max(1.0, maxAbsDiagonalFlat(workspace.matrix3, dimension)));

        double jitter = 0.0;
        double lowerBound = Double.NaN;
        for (int attempt = 0; attempt < 12; ++attempt) {
            GaussianMatrixOps.copyMatrixFlat(workspace.matrix3, workspace.matrix4, dimension);
            if (jitter > 0.0) {
                addDiagonalJitterFlat(workspace.matrix4, jitter, dimension);
            }
            logDet = invertPositiveDefiniteFromSymmetricCopyFlat(
                    workspace.matrix4, inverseOut, dimension, workspace.matrix5, workspace.matrix4);
            if (!Double.isNaN(logDet)) {
                return logDet;
            }
            if (jitter == 0.0) {
                if (Double.isNaN(lowerBound)) {
                    lowerBound = gershgorinLowerBoundFlat(workspace.matrix3, dimension);
                }
                jitter = lowerBound > 0.0 ? jitterBase : (-lowerBound + jitterBase);
            } else {
                jitter *= 10.0;
            }
        }
        throw new IllegalArgumentException("Matrix is not positive definite (failed robust inversion retries)");
    }

    private static double invertPositiveDefiniteFlatSubmatrix(final double[] matrix,
                                                               final double[] inverseOut,
                                                               final int dim,
                                                               final Workspace workspace) {
        System.arraycopy(matrix, 0, workspace.matrix3, 0, dim * dim);
        symmetrizeSquareFlat(workspace.matrix3, dim);

        double logDet = invertPositiveDefiniteFromSymmetricCopyFlat(
                workspace.matrix3, inverseOut, dim, workspace.matrix5, workspace.matrix4);
        if (!Double.isNaN(logDet)) {
            return logDet;
        }

        System.arraycopy(workspace.matrix3, 0, workspace.matrix4, 0, dim * dim);
        addDiagonalJitterFlat(workspace.matrix4, SYMMETRIC_JITTER_ABSOLUTE, dim);
        logDet = invertPositiveDefiniteFromSymmetricCopyFlat(
                workspace.matrix4, inverseOut, dim, workspace.matrix5, workspace.matrix4);
        if (!Double.isNaN(logDet)) {
            return logDet;
        }

        throw new IllegalArgumentException("Sub-matrix is not positive definite");
    }

    // Returns log-det(matrix) on success, Double.NaN if not positive definite.
    // Inlines Cholesky with simultaneous log-det accumulation to avoid a separate diagonal pass.
    private static double invertPositiveDefiniteFromSymmetricCopyFlat(final double[] matrix,
                                                                      final double[] inverseOut,
                                                                      final int dimension,
                                                                      final double[] choleskyScratch,
                                                                      final double[] lowerInverseScratch) {
        final int d = dimension;
        double logDet = 0.0;
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j <= i; ++j) {
                double sum = matrix[i * d + j];
                for (int k = 0; k < j; ++k) {
                    sum -= choleskyScratch[i * d + k] * choleskyScratch[j * d + k];
                }
                if (i == j) {
                    if (sum <= 0.0) return Double.NaN;
                    choleskyScratch[i * d + i] = Math.sqrt(sum);
                    logDet += Math.log(choleskyScratch[i * d + i]);
                } else {
                    final double denom = choleskyScratch[j * d + j];
                    if (denom == 0.0) return Double.NaN;
                    choleskyScratch[i * d + j] = sum / denom;
                }
            }
        }
        for (int col = 0; col < d; ++col) {
            for (int row = 0; row < d; ++row) {
                double sum = (row == col) ? 1.0 : 0.0;
                for (int k = 0; k < row; ++k) {
                    sum -= choleskyScratch[row * d + k] * lowerInverseScratch[k * d + col];
                }
                final double denom = choleskyScratch[row * d + row];
                if (denom == 0.0) return Double.NaN;
                lowerInverseScratch[row * d + col] = sum / denom;
            }
        }
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                double sum = 0.0;
                for (int k = 0; k < d; ++k) {
                    sum += lowerInverseScratch[k * d + i] * lowerInverseScratch[k * d + j];
                }
                inverseOut[i * d + j] = sum;
            }
        }
        GaussianMatrixOps.symmetrizeFlat(inverseOut, d);
        return 2.0 * logDet;
    }

    private static double observedQuadraticConstant(final CanonicalGaussianTransition transition,
                                                    final double[] observedValues,
                                                    final int[] observedIndices,
                                                    final int observedCount) {
        final int d = transition.getDimension();
        double quadratic = 0.0;
        double linear    = 0.0;
        for (int i = 0; i < observedCount; ++i) {
            final int obsI = observedIndices[i];
            final double valI = observedValues[obsI];
            linear += transition.informationY[obsI] * valI;
            for (int j = 0; j < observedCount; ++j) {
                final int obsJ = observedIndices[j];
                quadratic += valI * transition.precisionYY[obsI * d + obsJ] * observedValues[obsJ];
            }
        }
        return 0.5 * quadratic - linear;
    }

    private static double quadraticFormFlat(final double[] matrix, final double[] vector, final int dim) {
        double result = 0.0;
        for (int i = 0; i < dim; i++) {
            final int iOff = i * dim;
            double row = 0.0;
            for (int j = 0; j < dim; j++) {
                row += matrix[iOff + j] * vector[j];
            }
            result += vector[i] * row;
        }
        return result;
    }

    private static double dotFlat(final double[] left, final double[] right, final int dim) {
        double sum = 0.0;
        for (int i = 0; i < dim; i++) {
            sum += left[i] * right[i];
        }
        return sum;
    }

    private static void subtractMatrixInPlaceFlat(final double[] target, final double[] delta, final int dim) {
        final int d2 = dim * dim;
        for (int k = 0; k < d2; k++) {
            target[k] -= delta[k];
        }
    }

    // --- Block-extraction helpers (pair state has dimension 2d stored flat with stride 2d) ---

    private static void fillUpperLeftBlockFlat(final double[] source, final int d, final double[] out) {
        final int d2 = 2 * d;
        for (int i = 0; i < d; ++i) {
            System.arraycopy(source, i * d2, out, i * d, d);
        }
    }

    private static void fillUpperRightBlockFlat(final double[] source, final int d, final double[] out) {
        final int d2 = 2 * d;
        for (int i = 0; i < d; ++i) {
            System.arraycopy(source, i * d2 + d, out, i * d, d);
        }
    }

    private static void fillLowerLeftBlockFlat(final double[] source, final int d, final double[] out) {
        final int d2 = 2 * d;
        for (int i = 0; i < d; ++i) {
            System.arraycopy(source, (d + i) * d2, out, i * d, d);
        }
    }

    private static void fillLowerRightBlockFlat(final double[] source, final int d, final double[] out) {
        final int d2 = 2 * d;
        for (int i = 0; i < d; ++i) {
            System.arraycopy(source, (d + i) * d2 + d, out, i * d, d);
        }
    }

    private static void fillFirstBlockFlat(final double[] source, final int d, final double[] out) {
        System.arraycopy(source, 0, out, 0, d);
    }

    private static void fillSecondBlockFlat(final double[] source, final int d, final double[] out) {
        System.arraycopy(source, d, out, 0, d);
    }

    private static void ensureDimension(final Workspace workspace, final int dimension) {
        if (workspace.getDimension() != dimension) {
            throw new IllegalArgumentException(
                    "Workspace dimension mismatch: " + workspace.getDimension() + " vs " + dimension);
        }
    }

    private static void addDiagonalJitterFlat(final double[] matrix,
                                               final double jitter,
                                               final int dimension) {
        for (int i = 0; i < dimension; ++i) {
            matrix[i * dimension + i] += jitter;
        }
    }

    private static double maxAbsDiagonalFlat(final double[] matrix, final int dimension) {
        double max = 0.0;
        for (int i = 0; i < dimension; ++i) {
            max = Math.max(max, Math.abs(matrix[i * dimension + i]));
        }
        return max;
    }

    private static double gershgorinLowerBoundFlat(final double[] matrix, final int dimension) {
        double lowerBound = Double.POSITIVE_INFINITY;
        for (int i = 0; i < dimension; ++i) {
            final int iOff = i * dimension;
            double radius = 0.0;
            for (int j = 0; j < dimension; ++j) {
                if (i != j) radius += Math.abs(matrix[iOff + j]);
            }
            lowerBound = Math.min(lowerBound, matrix[iOff + i] - radius);
        }
        return lowerBound;
    }

    private static void symmetrizeSquareFlat(final double[] matrix, final int dimension) {
        for (int i = 0; i < dimension; ++i) {
            for (int j = i + 1; j < dimension; ++j) {
                final double avg = 0.5 * (matrix[i * dimension + j] + matrix[j * dimension + i]);
                matrix[i * dimension + j] = avg;
                matrix[j * dimension + i] = avg;
            }
        }
    }

    private static String formatVectorAsSquare(final double[] matrix, final int dimension) {
        final StringBuilder builder = new StringBuilder("[");
        for (int i = 0; i < dimension; ++i) {
            if (i > 0) builder.append(", ");
            builder.append('[');
            final int iOff = i * dimension;
            for (int j = 0; j < dimension; ++j) {
                if (j > 0) builder.append(", ");
                builder.append(matrix[iOff + j]);
            }
            builder.append(']');
        }
        builder.append(']');
        return builder.toString();
    }

    private static String formatVector(final double[] vector, final int dimension) {
        final StringBuilder builder = new StringBuilder("[");
        for (int i = 0; i < dimension; ++i) {
            if (i > 0) builder.append(", ");
            builder.append(vector[i]);
        }
        builder.append(']');
        return builder.toString();
    }
}
