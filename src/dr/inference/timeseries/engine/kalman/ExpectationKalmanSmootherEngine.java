package dr.inference.timeseries.engine.kalman;

import dr.evomodel.treedatalikelihood.continuous.canonical.math.MatrixOps;
import dr.inference.timeseries.core.TimeGrid;
import dr.inference.timeseries.model.gaussian.LinearGaussianObservationModel;
import dr.inference.timeseries.representation.GaussianTransitionRepresentation;

/**
 * Extension of the Kalman likelihood engine that appends an RTS backward smoothing pass.
 *
 * <p>After every forward pass the smoother computes, for every time step t = 0 … T−1:
 * <ul>
 *   <li>m_{t|T} = E[x_t | Y_{1:T}]          — smoothed posterior mean</li>
 *   <li>P_{t|T} = Cov(x_t | Y_{1:T})         — smoothed posterior covariance</li>
 *   <li>G_t    = P_{t|t} F_t^T P_{t+1|t}^{-1} — RTS smoother gain</li>
 * </ul>
 * The lag-1 cross-covariance P_{t+1,t|T} = P_{t+1|T} · G_t^T is <em>not</em> stored
 * explicitly; gradient formulas reconstruct it on demand from the stored gain and the
 * next step's smoothed covariance, saving O(T·d²) memory.
 *
 * <p>All results — including the forward trajectory — are cached behind a shared dirty
 * flag inherited from the base class. Calling {@link #makeDirty()} invalidates both the
 * likelihood and the smoother statistics atomically.
 */
public class ExpectationKalmanSmootherEngine extends ExpectationKalmanLikelihoodEngine implements GaussianSmootherResults {

    private final ForwardTrajectory trajectory;
    private final BranchSmootherStats[] smootherStats;

    // Pre-allocated working arrays for the backward pass — reused across calls.
    private final double[] backwardGain;
    private final double[] predictedCovCopy;
    private final double[] predictedCovCholesky;
    private final double[] predictedCovLowerInverse;
    private final double[] covDiff;
    private final double[]   meanDiff;
    private final double[] tempDxD1;
    private final double[] tempDxD2;
    private final double[]   tempD;

    public ExpectationKalmanSmootherEngine(final GaussianTransitionRepresentation transitionRepresentation,
                                final LinearGaussianObservationModel observationModel,
                                final TimeGrid timeGrid) {
        super(transitionRepresentation, observationModel, timeGrid);
        final int d = stateDimension;
        final int T = timeCount;

        trajectory    = new ForwardTrajectory(T, d);
        smootherStats = new BranchSmootherStats[T];
        for (int t = 0; t < T; ++t) {
            smootherStats[t] = new BranchSmootherStats(t, d, t < T - 1);
        }

        backwardGain      = new double[d * d];
        predictedCovCopy  = new double[d * d];
        predictedCovCholesky = new double[d * d];
        predictedCovLowerInverse = new double[d * d];
        covDiff           = new double[d * d];
        meanDiff          = new double[d];
        tempDxD1          = new double[d * d];
        tempDxD2          = new double[d * d];
        tempD             = new double[d];
    }

    // ─── Public accessors ───────────────────────────────────────────────────────

    /**
     * Returns the per-step smoother statistics, triggering a forward + backward pass
     * if the cache is stale. The returned array is owned by this engine — do not modify.
     */
    public BranchSmootherStats[] getSmootherStats() {
        getLogLikelihood();   // triggers computeLogLikelihood() if dirty
        return smootherStats;
    }

    /**
     * Returns the forward-pass trajectory snapshot, triggering a pass if stale.
     * The returned object is owned by this engine — do not modify.
     */
    public ForwardTrajectory getTrajectory() {
        getLogLikelihood();
        return trajectory;
    }

    public double[][] getSmoothedMeans() {
        getLogLikelihood();
        final double[][] out = new double[timeCount][stateDimension];
        for (int t = 0; t < timeCount; ++t) {
            ExpectationKalmanLikelihoodEngine.copyVector(smootherStats[t].smoothedMean, out[t]);
        }
        return out;
    }

    public double[][][] getSmoothedCovariances() {
        getLogLikelihood();
        final double[][][] out = new double[timeCount][stateDimension][stateDimension];
        for (int t = 0; t < timeCount; ++t) {
            MatrixOps.fromFlat(smootherStats[t].smoothedCovariance, out[t], stateDimension);
        }
        return out;
    }

    public double[][] getFilteredMeans() {
        getLogLikelihood();
        final double[][] out = new double[timeCount][stateDimension];
        for (int t = 0; t < timeCount; ++t) {
            trajectory.copyFilteredMeanTo(t, out[t]);
        }
        return out;
    }

    public double[][][] getFilteredCovariances() {
        getLogLikelihood();
        final double[][][] out = new double[timeCount][stateDimension][stateDimension];
        for (int t = 0; t < timeCount; ++t) {
            trajectory.copyFilteredCovarianceTo(t, out[t]);
        }
        return out;
    }

    public double[][] getPredictedMeans() {
        getLogLikelihood();
        final double[][] out = new double[timeCount][stateDimension];
        for (int t = 0; t < timeCount; ++t) {
            trajectory.copyPredictedMeanTo(t, out[t]);
        }
        return out;
    }

    public double[][][] getPredictedCovariances() {
        getLogLikelihood();
        final double[][][] out = new double[timeCount][stateDimension][stateDimension];
        for (int t = 0; t < timeCount; ++t) {
            trajectory.copyPredictedCovarianceTo(t, out[t]);
        }
        return out;
    }

    // ─── Override: run forward then backward ────────────────────────────────────

    /**
     * Runs the forward Kalman filter (storing full trajectory) followed by the RTS
     * backward smoother. Both results are cached together behind the parent's dirty flag.
     */
    @Override
    protected double computeLogLikelihood() {
        final double logLikelihood = runForwardPass(trajectory);
        runBackwardPass();
        return logLikelihood;
    }

    // ─── RTS backward pass ──────────────────────────────────────────────────────

    private void runBackwardPass() {
        final int T = trajectory.timeCount;
        final int d = trajectory.stateDimension;

        // Initialise: smoother at last step equals filter at last step
        System.arraycopy(
                trajectory.filteredMeans, trajectory.stateVectorOffset(T - 1),
                smootherStats[T - 1].smoothedMean, 0, d);
        System.arraycopy(
                trajectory.filteredCovariances, trajectory.stateMatrixOffset(T - 1),
                smootherStats[T - 1].smoothedCovariance, 0, trajectory.matrixSize);
        // smootherStats[T-1].smootherGain is null by construction — no forward step

        for (int t = T - 2; t >= 0; --t) {
            final int tVectorOffset = trajectory.stateVectorOffset(t);
            final int nextVectorOffset = trajectory.stateVectorOffset(t + 1);
            final int tMatrixOffset = trajectory.stateMatrixOffset(t);
            final int nextMatrixOffset = trajectory.stateMatrixOffset(t + 1);
            final int branchMatrixOffset = trajectory.branchMatrixOffset(t);

            // ── Smoother gain: G_t = P_{t|t} · F_t^T · P_{t+1|t}^{-1} ─────────
            //
            // We compute  P_{t|t} · F_t^T  first (using multiplyMatrixMatrixTransposedRight),
            // then right-multiply by P_{t+1|t}^{-1} obtained from a Cholesky factorisation.

            // Copy P_{t+1|t} before inverting (must not alias the stored trajectory).
            System.arraycopy(trajectory.predictedCovariances, nextMatrixOffset,
                    predictedCovCopy, 0, trajectory.matrixSize);
            if (!MatrixOps.tryCholesky(predictedCovCopy, predictedCovCholesky, d)) {
                throw new IllegalArgumentException("Matrix is not positive definite");
            }
            MatrixOps.invertFromCholesky(
                    predictedCovCholesky, predictedCovLowerInverse, predictedCovCopy, d);
            // predictedCovCopy now holds P_{t+1|t}^{-1}

            // tempDxD1 = P_{t|t} · F_t^T
            MatrixOps.matMulTransposedRight(
                    trajectory.filteredCovariances,
                    tMatrixOffset,
                    trajectory.transitionMatrices,
                    branchMatrixOffset,
                    tempDxD1,
                    d);

            // G_t = tempDxD1 · P_{t+1|t}^{-1}
            MatrixOps.matMul(tempDxD1, predictedCovCopy, backwardGain, d);
            System.arraycopy(backwardGain, 0, smootherStats[t].smootherGain, 0, trajectory.matrixSize);

            // ── Smoothed mean: m_{t|T} = m_{t|t} + G_t · (m_{t+1|T} − m_{t+1|t}) ──
            for (int i = 0; i < d; ++i) {
                meanDiff[i] = smootherStats[t + 1].smoothedMean[i]
                        - trajectory.predictedMeans[nextVectorOffset + i];
            }
            MatrixOps.matVec(backwardGain, meanDiff, tempD, d);
            for (int i = 0; i < d; ++i) {
                smootherStats[t].smoothedMean[i] =
                        trajectory.filteredMeans[tVectorOffset + i] + tempD[i];
            }

            // ── Smoothed covariance: P_{t|T} = P_{t|t} + G_t · (P_{t+1|T} − P_{t+1|t}) · G_t^T ──
            for (int i = 0; i < d; ++i) {
                final int rowOffset = i * d;
                for (int j = 0; j < d; ++j) {
                    covDiff[rowOffset + j] = smootherStats[t + 1].smoothedCovariance[rowOffset + j]
                            - trajectory.predictedCovariances[nextMatrixOffset + rowOffset + j];
                }
            }
            // tempDxD1 = G_t · (P_{t+1|T} − P_{t+1|t})
            MatrixOps.matMul(backwardGain, covDiff, tempDxD1, d);
            // tempDxD2 = tempDxD1 · G_t^T
            MatrixOps.matMulTransposedRight(tempDxD1, backwardGain, tempDxD2, d);
            // P_{t|T} = P_{t|t} + tempDxD2
            System.arraycopy(trajectory.filteredCovariances, tMatrixOffset,
                    smootherStats[t].smoothedCovariance, 0, trajectory.matrixSize);
            MatrixOps.addInPlace(smootherStats[t].smoothedCovariance, tempDxD2, trajectory.matrixSize);
            MatrixOps.symmetrize(smootherStats[t].smoothedCovariance, d);
        }
    }
}
