package dr.inference.timeseries.engine.kalman;

import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianMessageOps;

import dr.evomodel.treedatalikelihood.continuous.canonical.math.MatrixOps;
import dr.inference.timeseries.core.TimeGrid;
import dr.inference.timeseries.model.gaussian.LinearGaussianObservationModel;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianBranchTransitionKernel;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianState;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianTransition;
import dr.inference.timeseries.representation.CachedGaussianTransitionRepresentation;
import dr.inference.timeseries.representation.GaussianTransitionRepresentation;

/**
 * Canonical-form forward filter with RTS backward smoothing on the recovered moment trajectory.
 *
 * <p>The forward pass is performed in canonical form. At each step the predicted and
 * filtered canonical states are converted to means/covariances and stored in the same
 * trajectory structure used by the moment-form smoother. The backward pass then
 * applies the standard RTS recursions to those moments.
 */
public final class CanonicalKalmanSmootherEngine implements MomentSmootherResults, CanonicalSmootherResults {

    private final CanonicalGaussianBranchTransitionKernel canonicalKernel;
    private final GaussianTransitionRepresentation transitionRepresentation;
    private final CachedGaussianTransitionRepresentation cachedTransitionRepresentation;
    private final TimeGrid timeGrid;

    private final int stateDimension;
    private final int timeCount;

    private final CanonicalGaussianState futureMessage;
    private final CanonicalGaussianState backwardMessage;
    private final MomentForwardTrajectory trajectory;
    private final CanonicalForwardTrajectory canonicalTrajectory;
    private final CanonicalBranchGradientCache branchGradientCache;
    private final MomentBranchSmootherStats[] smootherStats;
    private final CanonicalGaussianMessageOps.Workspace messageWorkspace;
    private final CanonicalObservationUpdate observationUpdate;

    private final double[] stateVectorWorkspace;
    private final double[] flatMatrixWorkspace;
    private final double[] flatInverseWorkspace;
    private final double[] flatCholeskyWorkspace;
    private final double[] flatLowerInverseWorkspace;

    private boolean resultsKnown;
    private boolean momentTrajectoryKnown;
    private SharedCanonicalTimeSeriesSchedule sharedSchedule;
    private double logLikelihood;

    public CanonicalKalmanSmootherEngine(final CanonicalGaussianBranchTransitionKernel canonicalKernel,
                                         final GaussianTransitionRepresentation transitionRepresentation,
                                         final LinearGaussianObservationModel observationModel,
                                         final TimeGrid timeGrid) {
        if (canonicalKernel == null) {
            throw new IllegalArgumentException("canonicalKernel must not be null");
        }
        if (transitionRepresentation == null) {
            throw new IllegalArgumentException("transitionRepresentation must not be null");
        }
        if (observationModel == null) {
            throw new IllegalArgumentException("observationModel must not be null");
        }
        if (timeGrid == null) {
            throw new IllegalArgumentException("timeGrid must not be null");
        }
        this.canonicalKernel = canonicalKernel;
        this.transitionRepresentation = transitionRepresentation;
        this.cachedTransitionRepresentation = transitionRepresentation instanceof CachedGaussianTransitionRepresentation
                ? (CachedGaussianTransitionRepresentation) transitionRepresentation
                : null;
        this.timeGrid = timeGrid;
        this.stateDimension = canonicalKernel.getStateDimension();
        this.timeCount = timeGrid.getTimeCount();

        this.futureMessage = new CanonicalGaussianState(stateDimension);
        this.backwardMessage = new CanonicalGaussianState(stateDimension);
        this.trajectory = new MomentForwardTrajectory(timeCount, stateDimension);
        this.canonicalTrajectory = new CanonicalForwardTrajectory(timeCount, stateDimension);
        this.branchGradientCache = new CanonicalBranchGradientCache(
                timeCount, stateDimension, transitionRepresentation, timeGrid);
        this.smootherStats = new MomentBranchSmootherStats[timeCount];
        this.messageWorkspace = new CanonicalGaussianMessageOps.Workspace(stateDimension);
        for (int t = 0; t < timeCount; ++t) {
            smootherStats[t] = new MomentBranchSmootherStats(t, stateDimension, t < timeCount - 1);
        }

        this.observationUpdate = new CanonicalObservationUpdate(observationModel, stateDimension);
        this.stateVectorWorkspace = new double[stateDimension];
        this.flatMatrixWorkspace = new double[stateDimension * stateDimension];
        this.flatInverseWorkspace = new double[stateDimension * stateDimension];
        this.flatCholeskyWorkspace = new double[stateDimension * stateDimension];
        this.flatLowerInverseWorkspace = new double[stateDimension * stateDimension];
        this.sharedSchedule = null;
        if (cachedTransitionRepresentation != null) {
            cachedTransitionRepresentation.prepareTimeGrid(timeGrid);
        }
    }

    public double getLogLikelihood() {
        ensureResults();
        return logLikelihood;
    }

    @Override
    public void makeDirty() {
        resultsKnown = false;
        momentTrajectoryKnown = false;
        branchGradientCache.makeDirty();
        if (sharedSchedule != null) {
            sharedSchedule.makeDirty();
        }
    }

    public double[] getSmoothedMeansFlat() {
        ensureResults();
        ensureMomentTrajectory();
        final double[] out = new double[timeCount * stateDimension];
        for (int t = 0; t < timeCount; ++t) {
            System.arraycopy(smootherStats[t].smoothedMean, 0,
                    out, t * stateDimension, stateDimension);
        }
        return out;
    }

    public double[] getSmoothedCovariancesFlat() {
        ensureResults();
        ensureMomentTrajectory();
        final double[] out = new double[timeCount * stateDimension * stateDimension];
        final int matrixSize = stateDimension * stateDimension;
        for (int t = 0; t < timeCount; ++t) {
            System.arraycopy(smootherStats[t].smoothedCovariance, 0,
                    out, t * matrixSize, matrixSize);
        }
        return out;
    }

    public double[] getFilteredMeansFlat() {
        ensureResults();
        ensureMomentTrajectory();
        return trajectory.filteredMeans.clone();
    }

    public double[] getFilteredCovariancesFlat() {
        ensureResults();
        ensureMomentTrajectory();
        return trajectory.filteredCovariances.clone();
    }

    public double[] getPredictedMeansFlat() {
        ensureResults();
        ensureMomentTrajectory();
        return trajectory.predictedMeans.clone();
    }

    public double[] getPredictedCovariancesFlat() {
        ensureResults();
        ensureMomentTrajectory();
        return trajectory.predictedCovariances.clone();
    }

    private void ensureResults() {
        if (!resultsKnown) {
            logLikelihood = runForwardPass();
            runBackwardPass();
            resultsKnown = true;
            momentTrajectoryKnown = false;
        }
    }

    @Override
    public MomentBranchSmootherStats[] getSmootherStats() {
        ensureResults();
        ensureMomentTrajectory();
        return smootherStats;
    }

    @Override
    public MomentForwardTrajectory getTrajectory() {
        ensureResults();
        ensureMomentTrajectory();
        return trajectory;
    }

    @Override
    public CanonicalForwardTrajectory getCanonicalTrajectory() {
        ensureResults();
        return canonicalTrajectory;
    }

    @Override
    public CanonicalBranchGradientCache getCanonicalBranchGradientCache() {
        ensureResults();
        branchGradientCache.ensure(canonicalTrajectory);
        return branchGradientCache;
    }

    @Override
    public GaussianTransitionRepresentation getTransitionRepresentation() {
        return transitionRepresentation;
    }

    @Override
    public TimeGrid getTimeGrid() {
        return timeGrid;
    }

    @Override
    public int getTimeCount() {
        return timeCount;
    }

    @Override
    public int getStateDimension() {
        return stateDimension;
    }

    @Override
    public void setSharedSchedule(final SharedCanonicalTimeSeriesSchedule sharedSchedule) {
        this.sharedSchedule = sharedSchedule;
        branchGradientCache.setSharedSchedule(sharedSchedule);
    }

    private double runForwardPass() {
        observationUpdate.refreshStaticMatrices();

        double value = 0.0;
        CanonicalGaussianState previousFiltered = null;
        for (int timeIndex = 0; timeIndex < timeCount; ++timeIndex) {
            final CanonicalGaussianState predictedState =
                    canonicalTrajectory.predictedStates[timeIndex];
            final CanonicalGaussianState filteredState =
                    canonicalTrajectory.filteredStates[timeIndex];

            if (timeIndex == 0) {
                canonicalKernel.fillInitialCanonicalState(predictedState);
            } else {
                final CanonicalGaussianTransition transition =
                        canonicalTrajectory.transitions[timeIndex - 1];
                fillCanonicalTransition(timeIndex - 1, timeIndex, transition);
                predict(previousFiltered, transition, predictedState);
            }

            if (!observationUpdate.prepareForTime(timeIndex)) {
                copyState(predictedState, filteredState);
            } else {
                value += observationUpdate.applyTo(predictedState, filteredState, messageWorkspace);
            }
            previousFiltered = filteredState;
        }
        return value;
    }

    private void runBackwardPass() {
        final int T = trajectory.timeCount;
        copyState(canonicalTrajectory.filteredStates[T - 1], canonicalTrajectory.smoothedStates[T - 1]);

        for (int t = T - 2; t >= 0; --t) {
            subtractStates(canonicalTrajectory.smoothedStates[t + 1],
                    canonicalTrajectory.predictedStates[t + 1],
                    futureMessage);
            pushBackward(futureMessage,
                    canonicalTrajectory.transitions[t],
                    backwardMessage);
            combineStates(canonicalTrajectory.filteredStates[t],
                    backwardMessage,
                    canonicalTrajectory.smoothedStates[t]);
            buildPairPosterior(canonicalTrajectory.filteredStates[t],
                    canonicalTrajectory.transitions[t],
                    futureMessage,
                    canonicalTrajectory.branchPairStates[t]);
        }
    }

    private void ensureMomentTrajectory() {
        if (momentTrajectoryKnown) {
            return;
        }
        for (int timeIndex = 0; timeIndex < timeCount; ++timeIndex) {
            fillMomentsFromCanonical(canonicalTrajectory.predictedStates[timeIndex],
                    trajectory.predictedMeans,
                    trajectory.stateVectorOffset(timeIndex),
                    trajectory.predictedCovariances,
                    trajectory.stateMatrixOffset(timeIndex));
            fillMomentsFromCanonical(canonicalTrajectory.filteredStates[timeIndex],
                    trajectory.filteredMeans,
                    trajectory.stateVectorOffset(timeIndex),
                    trajectory.filteredCovariances,
                    trajectory.stateMatrixOffset(timeIndex));
            fillMomentsFromCanonical(canonicalTrajectory.smoothedStates[timeIndex],
                    smootherStats[timeIndex].smoothedMean,
                    smootherStats[timeIndex].smoothedCovariance);
            if (timeIndex < timeCount - 1) {
                transitionRepresentation.getTransitionMatrixFlat(timeIndex, timeIndex + 1, timeGrid,
                        flatMatrixWorkspace);
                trajectory.copyTransitionMatrixFrom(timeIndex, flatMatrixWorkspace);
                transitionRepresentation.getTransitionOffset(timeIndex, timeIndex + 1, timeGrid,
                        stateVectorWorkspace);
                trajectory.copyTransitionOffsetFrom(timeIndex, stateVectorWorkspace);
                transitionRepresentation.getTransitionCovarianceFlat(timeIndex, timeIndex + 1, timeGrid,
                        flatMatrixWorkspace);
                trajectory.copyStepCovarianceFrom(timeIndex, flatMatrixWorkspace);
            }
        }
        momentTrajectoryKnown = true;
    }

    private void predict(final CanonicalGaussianState previous,
                         final CanonicalGaussianTransition transition,
                         final CanonicalGaussianState out) {
        CanonicalGaussianMessageOps.pushForward(previous, transition, messageWorkspace, out);
    }

    private void fillMomentsFromCanonical(final CanonicalGaussianState canonical,
                                          final double[] meanOut,
                                          final double[] covarianceOut) {
        fillMomentsFromCanonical(canonical, meanOut, 0, covarianceOut, 0);
    }

    private void fillMomentsFromCanonical(final CanonicalGaussianState canonical,
                                          final double[] meanOut,
                                          final int meanOutOffset,
                                          final double[] covarianceOut,
                                          final int covarianceOutOffset) {
        invertPositiveDefiniteFlatInput(canonical.precision, covarianceOut, covarianceOutOffset, stateDimension);
        MatrixOps.matVec(covarianceOut, covarianceOutOffset, canonical.information,
                stateVectorWorkspace, stateDimension);
        System.arraycopy(stateVectorWorkspace, 0, meanOut, meanOutOffset, stateDimension);
    }

    private void invertPositiveDefiniteFlatInput(final double[] matrix,
                                                final double[] inverseOut,
                                                final int inverseOutOffset,
                                                final int dimension) {
        if (!MatrixOps.tryCholesky(matrix, flatCholeskyWorkspace, dimension)) {
            throw new IllegalArgumentException("Matrix is not positive definite");
        }
        MatrixOps.invertFromCholesky(flatCholeskyWorkspace, flatLowerInverseWorkspace, flatInverseWorkspace, dimension);
        System.arraycopy(flatInverseWorkspace, 0, inverseOut, inverseOutOffset, dimension * dimension);
    }

    private void fillCanonicalTransition(final int fromIndex,
                                         final int toIndex,
                                         final CanonicalGaussianTransition out) {
        if (cachedTransitionRepresentation != null) {
            cachedTransitionRepresentation.getCanonicalTransition(fromIndex, toIndex, timeGrid, out);
        } else {
            canonicalKernel.fillCanonicalTransition(validatedDelta(fromIndex, toIndex), out);
        }
    }

    private double validatedDelta(final int fromIndex, final int toIndex) {
        final double dt = timeGrid.getDelta(fromIndex, toIndex);
        if (!(dt > 0.0)) {
            throw new IllegalArgumentException("Time increments must be strictly positive");
        }
        return dt;
    }

    private static void copyState(final CanonicalGaussianState source, final CanonicalGaussianState target) {
        System.arraycopy(source.precision, 0, target.precision, 0, source.precision.length);
        System.arraycopy(source.information, 0, target.information, 0, source.information.length);
        target.logNormalizer = source.logNormalizer;
    }

    private void subtractStates(final CanonicalGaussianState left,
                                final CanonicalGaussianState right,
                                final CanonicalGaussianState out) {
        for (int i = 0; i < stateDimension; ++i) {
            out.information[i] = left.information[i] - right.information[i];
            for (int j = 0; j < stateDimension; ++j) {
                out.precision[i * stateDimension + j] = left.precision[i * stateDimension + j] - right.precision[i * stateDimension + j];
            }
        }
        out.logNormalizer = left.logNormalizer - right.logNormalizer;
    }

    private void buildPairPosterior(final CanonicalGaussianState filteredState,
                                    final CanonicalGaussianTransition transition,
                                    final CanonicalGaussianState futureMessage,
                                    final CanonicalGaussianState pairOut) {
        CanonicalGaussianMessageOps.buildPairPosterior(filteredState, transition, futureMessage, pairOut);
    }

    private void pushBackward(final CanonicalGaussianState futureMessage,
                              final CanonicalGaussianTransition transition,
                              final CanonicalGaussianState out) {
        CanonicalGaussianMessageOps.pushBackward(futureMessage, transition, messageWorkspace, out);
    }

    private static void combineStates(final CanonicalGaussianState left,
                                      final CanonicalGaussianState right,
                                      final CanonicalGaussianState out) {
        CanonicalGaussianMessageOps.combineStates(left, right, out);
    }

}
