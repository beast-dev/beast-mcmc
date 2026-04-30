package dr.inference.timeseries.engine.gaussian;

import dr.evomodel.treedatalikelihood.continuous.gaussian.message.GaussianMatrixOps;

import dr.inference.timeseries.core.TimeGrid;
import dr.inference.timeseries.engine.LikelihoodEngine;
import dr.inference.timeseries.gaussian.GaussianObservationModel;
import dr.inference.timeseries.representation.GaussianTransitionRepresentation;

/**
 * Linear-Gaussian likelihood engine based on a Kalman filter.
 *
 * <p>This engine marginalizes over the latent state trajectory. Missing observations
 * are supported by encoding an entire observation column as {@link Double#NaN}.
 *
 * <p>The forward pass is factored into {@link #runForwardPass(ForwardTrajectory)} so
 * that subclasses (e.g. {@link KalmanSmootherEngine}) can store the full trajectory
 * for a subsequent RTS backward pass without duplicating the filter logic. When
 * {@code trajectoryOut} is {@code null} the trajectory is computed but not stored,
 * keeping the base-class path allocation-free.
 *
 * <p>All static linear-algebra helpers are package-private so that other classes in
 * this package (smoother, gradient formulas) can reuse them without duplication.
 */
public class KalmanLikelihoodEngine implements LikelihoodEngine {

    static final double LOG_TWO_PI = GaussianMatrixOps.LOG_TWO_PI;
    static final double MIN_DIAGONAL_JITTER = GaussianMatrixOps.MIN_DIAGONAL_JITTER;

    private final GaussianTransitionRepresentation transitionRepresentation;
    private final GaussianObservationModel observationModel;
    private final TimeGrid timeGrid;

    final int stateDimension;
    final int observationDimension;
    final int timeCount;

    private final double[] filteredMean;
    private final double[] predictedMean;
    private final double[] offset;
    private final double[] innovation;
    private final double[] kalmanGainColumn;
    private final double[] workingVector;
    private final double[] observationVector;

    private final double[][] filteredCovariance;
    private final double[][] predictedCovariance;
    private final double[][] transitionMatrix;
    private final double[][] transitionCovariance;
    private final double[][] initialCovariance;
    private final double[][] designMatrix;
    private final double[][] noiseCovariance;
    private final double[][] innovationCovariance;
    private final double[][] innovationCovarianceInverse;
    private final double[][] gainMatrix;
    private final double[][] tempStateState;
    private final double[][] tempStateObs;
    private final double[][] tempObsState;
    private final double[][] tempObsObs;

    private boolean likelihoodKnown;
    private double logLikelihood;

    public KalmanLikelihoodEngine(final GaussianTransitionRepresentation transitionRepresentation,
                                  final GaussianObservationModel observationModel,
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
        kalmanGainColumn          = new double[stateDimension];
        workingVector             = new double[stateDimension];
        observationVector         = new double[observationDimension];

        filteredCovariance        = new double[stateDimension][stateDimension];
        predictedCovariance       = new double[stateDimension][stateDimension];
        transitionMatrix          = new double[stateDimension][stateDimension];
        transitionCovariance      = new double[stateDimension][stateDimension];
        initialCovariance         = new double[stateDimension][stateDimension];
        designMatrix              = new double[observationDimension][stateDimension];
        noiseCovariance           = new double[observationDimension][observationDimension];
        innovationCovariance      = new double[observationDimension][observationDimension];
        innovationCovarianceInverse = new double[observationDimension][observationDimension];
        gainMatrix                = new double[stateDimension][observationDimension];
        tempStateState            = new double[stateDimension][stateDimension];
        tempStateObs              = new double[stateDimension][observationDimension];
        tempObsState              = new double[observationDimension][stateDimension];
        tempObsObs                = new double[observationDimension][observationDimension];
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

    public GaussianObservationModel getObservationModel() {
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
        observationModel.fillDesignMatrix(designMatrix);
        observationModel.fillNoiseCovariance(noiseCovariance);
        transitionRepresentation.getInitialMean(filteredMean);
        transitionRepresentation.getInitialCovariance(initialCovariance);
        copyMatrix(initialCovariance, filteredCovariance);

        double value = 0.0;

        for (int timeIndex = 0; timeIndex < timeCount; ++timeIndex) {

            // ── Prediction ──────────────────────────────────────────────────────
            if (timeIndex == 0) {
                copyVector(filteredMean, predictedMean);
                copyMatrix(filteredCovariance, predictedCovariance);
            } else {
                transitionRepresentation.getTransitionMatrix(
                        timeIndex - 1, timeIndex, timeGrid, transitionMatrix);
                transitionRepresentation.getTransitionOffset(
                        timeIndex - 1, timeIndex, timeGrid, offset);
                transitionRepresentation.getTransitionCovariance(
                        timeIndex - 1, timeIndex, timeGrid, transitionCovariance);

                multiplyMatrixVector(transitionMatrix, filteredMean, predictedMean);
                addInPlace(predictedMean, offset);

                multiplyMatrixMatrix(transitionMatrix, filteredCovariance, tempStateState);
                multiplyMatrixMatrixTransposedRight(tempStateState, transitionMatrix, predictedCovariance);
                addMatrixInPlace(predictedCovariance, transitionCovariance);
                symmetrize(predictedCovariance);

                if (trajectoryOut != null) {
                    copyMatrix(transitionMatrix,    trajectoryOut.transitionMatrices[timeIndex - 1]);
                    copyVector(offset,              trajectoryOut.transitionOffsets[timeIndex - 1]);
                    copyMatrix(transitionCovariance, trajectoryOut.stepCovariances[timeIndex - 1]);
                }
            }

            if (trajectoryOut != null) {
                copyVector(predictedMean,       trajectoryOut.predictedMeans[timeIndex]);
                copyMatrix(predictedCovariance, trajectoryOut.predictedCovariances[timeIndex]);
            }

            // ── Update ──────────────────────────────────────────────────────────
            if (observationModel.isObservationMissing(timeIndex)) {
                copyVector(predictedMean,       filteredMean);
                copyMatrix(predictedCovariance, filteredCovariance);
            } else {
                observationModel.fillObservationVector(timeIndex, observationVector);

                multiplyMatrixVector(designMatrix, predictedMean, innovation);
                for (int i = 0; i < observationDimension; ++i) {
                    innovation[i] = observationVector[i] - innovation[i];
                }

                multiplyMatrixMatrix(designMatrix, predictedCovariance, tempObsState);
                multiplyMatrixMatrixTransposedRight(tempObsState, designMatrix, innovationCovariance);
                addMatrixInPlace(innovationCovariance, noiseCovariance);
                symmetrize(innovationCovariance);
                addJitterToDiagonal(innovationCovariance, MIN_DIAGONAL_JITTER);

                final CholeskyFactor innovationCholesky = cholesky(innovationCovariance);
                final double logDetInnovation = innovationCholesky.logDeterminant();
                copyMatrix(innovationCovariance, innovationCovarianceInverse);
                invertPositiveDefiniteFromCholesky(innovationCovarianceInverse, innovationCholesky);

                value += -0.5 * (observationDimension * LOG_TWO_PI
                        + logDetInnovation
                        + quadraticForm(innovationCovarianceInverse, innovation));

                multiplyMatrixMatrixTransposedRight(predictedCovariance, designMatrix, tempStateObs);
                multiplyMatrixMatrix(tempStateObs, innovationCovarianceInverse, gainMatrix,
                        stateDimension, observationDimension, observationDimension);

                multiplyMatrixVector(gainMatrix, innovation, workingVector, stateDimension, observationDimension);
                for (int i = 0; i < stateDimension; ++i) {
                    filteredMean[i] = predictedMean[i] + workingVector[i];
                }

                // Joseph-stabilised covariance update: P = (I−KH) P (I−KH)^T + K R K^T
                //
                // predictedCovariance is no longer needed after the gain computation above,
                // so it is reused as a [d×d] scratch buffer. This avoids the shape mismatch
                // that would arise from using tempObsObs ([p×p]) for a [d×d] intermediate.
                multiplyMatrixMatrix(gainMatrix, designMatrix, tempStateState,
                        stateDimension, observationDimension, stateDimension);
                identityMinus(tempStateState);
                multiplyMatrixMatrix(tempStateState, predictedCovariance, filteredCovariance,
                        stateDimension, stateDimension, stateDimension);
                multiplyMatrixMatrixTransposedRight(filteredCovariance, tempStateState, predictedCovariance);
                copyMatrix(predictedCovariance, filteredCovariance);
                multiplyMatrixMatrix(gainMatrix, noiseCovariance, tempStateObs,
                        stateDimension, observationDimension, observationDimension);
                multiplyMatrixMatrixTransposedRight(tempStateObs, gainMatrix, tempStateState);
                addMatrixInPlace(filteredCovariance, tempStateState);
                symmetrize(filteredCovariance);
            }

            if (trajectoryOut != null) {
                copyVector(filteredMean,       trajectoryOut.filteredMeans[timeIndex]);
                copyMatrix(filteredCovariance, trajectoryOut.filteredCovariances[timeIndex]);
            }
        }

        return value;
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
