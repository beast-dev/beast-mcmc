package dr.inference.timeseries.engine.kalman;

import dr.evomodel.treedatalikelihood.continuous.canonical.math.MatrixOps;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.GaussianMatrixOps;

import dr.inference.timeseries.core.TimeGrid;
import dr.inference.timeseries.engine.LikelihoodEngine;
import dr.inference.timeseries.model.gaussian.LinearGaussianObservationModel;
import dr.inference.timeseries.representation.GaussianTransitionRepresentation;

/**
 * Linear-Gaussian likelihood engine based on a Kalman filter.
 *
 * <p>This engine marginalizes over the latent state trajectory. Missing observations
 * are supported by encoding an entire observation column as {@link Double#NaN}.
 *
 * <p>The forward pass is factored into {@link #runForwardPass(ForwardTrajectory)} so
 * that subclasses (e.g. {@link ExpectationKalmanSmootherEngine}) can store the full trajectory
 * for a subsequent RTS backward pass without duplicating the filter logic. When
 * {@code trajectoryOut} is {@code null} the trajectory is computed but not stored,
 * keeping the base-class path allocation-free.
 *
 * <p>All static linear-algebra helpers are package-private so that other classes in
 * this package (smoother, gradient formulas) can reuse them without duplication.
 */
public class ExpectationKalmanLikelihoodEngine implements LikelihoodEngine {

    static final double LOG_TWO_PI = GaussianMatrixOps.LOG_TWO_PI;
    static final double MIN_DIAGONAL_JITTER = GaussianMatrixOps.MIN_DIAGONAL_JITTER;

    private final GaussianTransitionRepresentation transitionRepresentation;
    private final LinearGaussianObservationModel observationModel;
    private final TimeGrid timeGrid;

    final int stateDimension;
    final int observationDimension;
    final int timeCount;

    private final double[] filteredMean;
    private final double[] predictedMean;
    private final double[] offset;
    private final double[] innovation;
    private final double[] workingVector;
    private final double[] observationVector;
    private final double[] observationWorkingVector;

    private final double[] filteredCovariance;
    private final double[] predictedCovariance;
    private final double[] transitionMatrix;
    private final double[] transitionCovariance;
    private final double[] initialCovariance;
    private final double[] designMatrix;
    private final double[] noiseCovariance;
    private final double[] innovationCovariance;
    private final double[] innovationCovarianceInverse;
    private final double[] innovationCholeskyLower;
    private final double[] innovationLowerInverse;
    private final double[] gainMatrix;
    private final double[] tempStateState;
    private final double[] tempStateObs;
    private final double[] tempObsState;

    private boolean likelihoodKnown;
    private double logLikelihood;

    public ExpectationKalmanLikelihoodEngine(final GaussianTransitionRepresentation transitionRepresentation,
                                  final LinearGaussianObservationModel observationModel,
                                  final TimeGrid timeGrid) {
        if (transitionRepresentation == null) {
            throw new IllegalArgumentException("transitionRepresentation must not be null");
        }
        if (observationModel == null) {
            throw new IllegalArgumentException("observationModel must not be null");
        }
        if (timeGrid == null) {
            throw new IllegalArgumentException("timeGrid must not be null");
        }
        this.transitionRepresentation = transitionRepresentation;
        this.observationModel         = observationModel;
        this.timeGrid                 = timeGrid;
        this.stateDimension           = transitionRepresentation.getStateDimension();
        this.observationDimension     = observationModel.getObservationDimension();
        this.timeCount                = timeGrid.getTimeCount();

        if (observationModel.getTimeCount() != timeCount) {
            throw new IllegalArgumentException(
                    "Observation time count (" + observationModel.getTimeCount() +
                    ") must match time grid count (" + timeCount + ")");
        }

        filteredMean              = new double[stateDimension];
        predictedMean             = new double[stateDimension];
        offset                    = new double[stateDimension];
        innovation                = new double[observationDimension];
        workingVector             = new double[stateDimension];
        observationVector         = new double[observationDimension];
        observationWorkingVector  = new double[observationDimension];

        filteredCovariance        = new double[stateDimension * stateDimension];
        predictedCovariance       = new double[stateDimension * stateDimension];
        transitionMatrix          = new double[stateDimension * stateDimension];
        transitionCovariance      = new double[stateDimension * stateDimension];
        initialCovariance         = new double[stateDimension * stateDimension];
        designMatrix              = new double[observationDimension * stateDimension];
        noiseCovariance           = new double[observationDimension * observationDimension];
        innovationCovariance      = new double[observationDimension * observationDimension];
        innovationCovarianceInverse = new double[observationDimension * observationDimension];
        innovationCholeskyLower   = new double[observationDimension * observationDimension];
        innovationLowerInverse    = new double[observationDimension * observationDimension];
        gainMatrix                = new double[stateDimension * observationDimension];
        tempStateState            = new double[stateDimension * stateDimension];
        tempStateObs              = new double[stateDimension * observationDimension];
        tempObsState              = new double[observationDimension * stateDimension];
    }

    // ─── LikelihoodEngine ───────────────────────────────────────────────────────

    @Override
    public double getLogLikelihood() {
        if (!likelihoodKnown) {
            logLikelihood = computeLogLikelihood();
            likelihoodKnown = true;
        }
        return logLikelihood;
    }

    @Override
    public void makeDirty() {
        likelihoodKnown = false;
    }

    // ─── Accessors ──────────────────────────────────────────────────────────────

    public GaussianTransitionRepresentation getTransitionRepresentation() {
        return transitionRepresentation;
    }

    public LinearGaussianObservationModel getObservationModel() {
        return observationModel;
    }

    public TimeGrid getTimeGrid() {
        return timeGrid;
    }

    // ─── Forward pass ───────────────────────────────────────────────────────────

    /**
     * Hook for subclasses: by default delegates to {@link #runForwardPass(ForwardTrajectory)}
     * with a null trajectory (no storage overhead).
     */
    protected double computeLogLikelihood() {
        return runForwardPass(null);
    }

    /**
     * Runs the Kalman filter forward pass and returns the log-likelihood.
     *
     * <p>When {@code trajectoryOut} is non-null, filtered/predicted moments and transition
     * parameters are written into it at every step. This is the only place where trajectory
     * storage is performed; the rest of the engine is unaffected.
     *
     * @param trajectoryOut destination for the full trajectory, or {@code null} to skip storage
     * @return log p(Y | θ) accumulated over all observed steps
     */
    protected double runForwardPass(final ForwardTrajectory trajectoryOut) {
        observationModel.fillDesignMatrixFlat(designMatrix, stateDimension);
        observationModel.fillNoiseCovarianceFlat(noiseCovariance);
        transitionRepresentation.getInitialMean(filteredMean);
        transitionRepresentation.getInitialCovarianceFlat(initialCovariance);
        copyFlatMatrix(initialCovariance, filteredCovariance, stateDimension);

        double value = 0.0;

        for (int timeIndex = 0; timeIndex < timeCount; ++timeIndex) {

            // ── Prediction ──────────────────────────────────────────────────────
            if (timeIndex == 0) {
                copyVector(filteredMean, predictedMean);
                copyFlatMatrix(filteredCovariance, predictedCovariance, stateDimension);
            } else {
                transitionRepresentation.getTransitionMatrixFlat(
                        timeIndex - 1, timeIndex, timeGrid, transitionMatrix);
                transitionRepresentation.getTransitionOffset(
                        timeIndex - 1, timeIndex, timeGrid, offset);
                transitionRepresentation.getTransitionCovarianceFlat(
                        timeIndex - 1, timeIndex, timeGrid, transitionCovariance);

                MatrixOps.matVec(transitionMatrix, filteredMean, predictedMean, stateDimension);
                addInPlace(predictedMean, offset);

                MatrixOps.matMul(transitionMatrix, filteredCovariance, tempStateState, stateDimension);
                MatrixOps.matMulTransposedRight(tempStateState, transitionMatrix, predictedCovariance, stateDimension);
                MatrixOps.addInPlace(predictedCovariance, transitionCovariance, stateDimension * stateDimension);
                MatrixOps.symmetrize(predictedCovariance, stateDimension);

                if (trajectoryOut != null) {
                    MatrixOps.fromFlat(transitionMatrix, trajectoryOut.transitionMatrices[timeIndex - 1], stateDimension);
                    copyVector(offset,              trajectoryOut.transitionOffsets[timeIndex - 1]);
                    MatrixOps.fromFlat(transitionCovariance, trajectoryOut.stepCovariances[timeIndex - 1], stateDimension);
                }
            }

            if (trajectoryOut != null) {
                copyVector(predictedMean,       trajectoryOut.predictedMeans[timeIndex]);
                MatrixOps.fromFlat(predictedCovariance, trajectoryOut.predictedCovariances[timeIndex], stateDimension);
            }

            // ── Update ──────────────────────────────────────────────────────────
            if (observationModel.isObservationMissing(timeIndex)) {
                copyVector(predictedMean,       filteredMean);
                copyFlatMatrix(predictedCovariance, filteredCovariance, stateDimension);
            } else {
                observationModel.fillObservationVector(timeIndex, observationVector);

                MatrixOps.matVec(designMatrix, predictedMean, innovation, observationDimension, stateDimension);
                for (int i = 0; i < observationDimension; ++i) {
                    innovation[i] = observationVector[i] - innovation[i];
                }

                MatrixOps.matMul(designMatrix, predictedCovariance, tempObsState,
                        observationDimension, stateDimension, stateDimension);
                multiplyFlatTransposedRight(
                        tempObsState, designMatrix, innovationCovariance,
                        observationDimension, stateDimension, observationDimension);
                MatrixOps.addInPlace(innovationCovariance, noiseCovariance,
                        observationDimension * observationDimension);
                MatrixOps.symmetrize(innovationCovariance, observationDimension);
                addJitterToDiagonalFlat(innovationCovariance, observationDimension, MIN_DIAGONAL_JITTER);

                if (!MatrixOps.tryCholesky(
                        innovationCovariance, innovationCholeskyLower, observationDimension)) {
                    throw new IllegalArgumentException("Matrix is not positive definite");
                }
                final double logDetInnovation = MatrixOps.invertFromCholesky(
                        innovationCholeskyLower,
                        innovationLowerInverse,
                        innovationCovarianceInverse,
                        observationDimension);

                value += -0.5 * (observationDimension * LOG_TWO_PI
                        + logDetInnovation
                        + MatrixOps.quadraticForm(
                        innovation, innovationCovarianceInverse, observationDimension, observationWorkingVector));

                multiplyFlatTransposedRight(
                        predictedCovariance, designMatrix, tempStateObs,
                        stateDimension, stateDimension, observationDimension);
                MatrixOps.matMul(tempStateObs, innovationCovarianceInverse, gainMatrix,
                        stateDimension, observationDimension, observationDimension);

                MatrixOps.matVec(gainMatrix, innovation, workingVector, stateDimension, observationDimension);
                for (int i = 0; i < stateDimension; ++i) {
                    filteredMean[i] = predictedMean[i] + workingVector[i];
                }

                // Joseph-stabilised covariance update: P = (I−KH) P (I−KH)^T + K R K^T
                //
                // predictedCovariance is no longer needed after the gain computation above,
                // so it is reused as a [d x d] scratch buffer.
                MatrixOps.matMul(gainMatrix, designMatrix, tempStateState,
                        stateDimension, observationDimension, stateDimension);
                identityMinusFlat(tempStateState, stateDimension);
                MatrixOps.matMul(tempStateState, predictedCovariance, filteredCovariance,
                        stateDimension, stateDimension, stateDimension);
                MatrixOps.matMulTransposedRight(filteredCovariance, tempStateState, predictedCovariance, stateDimension);
                copyFlatMatrix(predictedCovariance, filteredCovariance, stateDimension);
                MatrixOps.matMul(gainMatrix, noiseCovariance, tempStateObs,
                        stateDimension, observationDimension, observationDimension);
                multiplyFlatTransposedRight(
                        tempStateObs, gainMatrix, tempStateState,
                        stateDimension, observationDimension, stateDimension);
                MatrixOps.addInPlace(filteredCovariance, tempStateState, stateDimension * stateDimension);
                MatrixOps.symmetrize(filteredCovariance, stateDimension);
            }

            if (trajectoryOut != null) {
                copyVector(filteredMean,       trajectoryOut.filteredMeans[timeIndex]);
                MatrixOps.fromFlat(filteredCovariance, trajectoryOut.filteredCovariances[timeIndex], stateDimension);
            }
        }

        return value;
    }

    private static void copyFlatMatrix(final double[] source,
                                       final double[] target,
                                       final int dimension) {
        System.arraycopy(source, 0, target, 0, dimension * dimension);
    }

    /**
     * Computes {@code out = left * right^T} for row-major rectangular matrices.
     * left is rows x inner, right is cols x inner, and out is rows x cols.
     */
    private static void multiplyFlatTransposedRight(final double[] left,
                                                   final double[] right,
                                                   final double[] out,
                                                   final int rows,
                                                   final int inner,
                                                   final int cols) {
        for (int i = 0; i < rows; ++i) {
            final int leftRowOffset = i * inner;
            final int outRowOffset = i * cols;
            for (int j = 0; j < cols; ++j) {
                final int rightRowOffset = j * inner;
                double sum = 0.0;
                for (int k = 0; k < inner; ++k) {
                    sum += left[leftRowOffset + k] * right[rightRowOffset + k];
                }
                out[outRowOffset + j] = sum;
            }
        }
    }

    private static void identityMinusFlat(final double[] matrix,
                                          final int dimension) {
        for (int i = 0; i < dimension; ++i) {
            final int rowOffset = i * dimension;
            for (int j = 0; j < dimension; ++j) {
                matrix[rowOffset + j] = (i == j ? 1.0 : 0.0) - matrix[rowOffset + j];
            }
        }
    }

    private static void addJitterToDiagonalFlat(final double[] matrix,
                                                final int dimension,
                                                final double minimumJitter) {
        for (int i = 0; i < dimension; ++i) {
            final int index = i * dimension + i;
            if (matrix[index] < minimumJitter) {
                matrix[index] = minimumJitter;
            }
        }
    }

    // ─── Package-private linear-algebra compatibility wrappers ───────────────────
    // These stay here so existing expectation-form code in this package can keep
    // compiling while the actual implementations live in GaussianMatrixOps.

    static void multiplyMatrixVector(final double[][] matrix,
                                     final double[] vector,
                                     final double[] out) {
        GaussianMatrixOps.multiplyMatrixVector(matrix, vector, out);
    }

    static void multiplyMatrixVector(final double[][] matrix,
                                     final double[] vector,
                                     final double[] out,
                                     final int rows,
                                     final int cols) {
        GaussianMatrixOps.multiplyMatrixVector(matrix, vector, out, rows, cols);
    }

    static void multiplyMatrixMatrix(final double[][] left,
                                     final double[][] right,
                                     final double[][] out) {
        GaussianMatrixOps.multiplyMatrixMatrix(left, right, out);
    }

    static void multiplyMatrixMatrix(final double[][] left,
                                     final double[][] right,
                                     final double[][] out,
                                     final int rows,
                                     final int inner,
                                     final int cols) {
        GaussianMatrixOps.multiplyMatrixMatrix(left, right, out, rows, inner, cols);
    }

    /** Computes {@code out = left · right^T}. */
    static void multiplyMatrixMatrixTransposedRight(final double[][] left,
                                                    final double[][] right,
                                                    final double[][] out) {
        GaussianMatrixOps.multiplyMatrixMatrixTransposedRight(left, right, out);
    }

    static void addMatrixInPlace(final double[][] accumulator, final double[][] increment) {
        GaussianMatrixOps.addMatrixInPlace(accumulator, increment);
    }

    static void addInPlace(final double[] accumulator, final double[] increment) {
        GaussianMatrixOps.addInPlace(accumulator, increment);
    }

    static void identityMinus(final double[][] matrix) {
        GaussianMatrixOps.identityMinus(matrix);
    }

    static void copyVector(final double[] source, final double[] target) {
        GaussianMatrixOps.copyVector(source, target);
    }

    static void copyMatrix(final double[][] source, final double[][] target) {
        GaussianMatrixOps.copyMatrix(source, target);
    }

    static void copyMatrix(final double[][] source,
                           final double[][] target,
                           final int rows,
                           final int cols) {
        GaussianMatrixOps.copyMatrix(source, target, rows, cols);
    }

    static void symmetrize(final double[][] matrix) {
        GaussianMatrixOps.symmetrize(matrix);
    }

    static void addJitterToDiagonal(final double[][] matrix, final double minimumJitter) {
        GaussianMatrixOps.addJitterToDiagonal(matrix, minimumJitter);
    }

    static double quadraticForm(final double[][] matrix, final double[] vector) {
        return GaussianMatrixOps.quadraticForm(matrix, vector);
    }

    static CholeskyFactor cholesky(final double[][] matrix) {
        return new CholeskyFactor(GaussianMatrixOps.cholesky(matrix));
    }

    static void invertPositiveDefiniteFromCholesky(final double[][] out,
                                                   final CholeskyFactor factor) {
        GaussianMatrixOps.invertPositiveDefiniteFromCholesky(out, factor.delegate);
    }

    // ─── Cholesky factor ────────────────────────────────────────────────────────

    static final class CholeskyFactor {
        private final GaussianMatrixOps.CholeskyFactor delegate;

        CholeskyFactor(final GaussianMatrixOps.CholeskyFactor delegate) {
            this.delegate = delegate;
        }

        double logDeterminant() {
            return delegate.logDeterminant();
        }

        void solveSymmetricSystem(final double[] rhs, final double[] out) {
            delegate.solveSymmetricSystem(rhs, out);
        }
    }
}
