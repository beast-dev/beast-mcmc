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
 * trajectory structure used by the expectation-form smoother. The backward pass then
 * applies the standard RTS recursions to those moments.
 */
public final class CanonicalKalmanSmootherEngine implements GaussianSmootherResults, CanonicalSmootherResults {

    private final CanonicalGaussianBranchTransitionKernel canonicalKernel;
    private final GaussianTransitionRepresentation transitionRepresentation;
    private final CachedGaussianTransitionRepresentation cachedTransitionRepresentation;
    private final LinearGaussianObservationModel observationModel;
    private final TimeGrid timeGrid;

    private final int stateDimension;
    private final int timeCount;
    private final int observationDimension;
    private final int maximumMatrixDimension;

    private final CanonicalGaussianState futureMessage;
    private final CanonicalGaussianState backwardMessage;
    private final ForwardTrajectory trajectory;
    private final CanonicalForwardTrajectory canonicalTrajectory;
    private final CanonicalBranchGradientCache branchGradientCache;
    private final BranchSmootherStats[] smootherStats;
    private final CanonicalGaussianMessageOps.Workspace messageWorkspace;

    private final double[] designMatrix;
    private final double[] noiseCovariance;
    private final double[] noisePrecision;
    private final double[] observationPrecisionContribution;
    private final double[] obsWorkspace;
    private final double[] observationVector;
    private final double[] observationInformation;
    private final double[] observationVectorWorkspace;
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
        this.observationModel = observationModel;
        this.timeGrid = timeGrid;
        this.stateDimension = canonicalKernel.getStateDimension();
        this.timeCount = timeGrid.getTimeCount();
        this.observationDimension = observationModel.getObservationDimension();
        this.maximumMatrixDimension = Math.max(stateDimension, observationDimension);

        this.futureMessage = new CanonicalGaussianState(stateDimension);
        this.backwardMessage = new CanonicalGaussianState(stateDimension);
        this.trajectory = new ForwardTrajectory(timeCount, stateDimension);
        this.canonicalTrajectory = new CanonicalForwardTrajectory(timeCount, stateDimension);
        this.branchGradientCache = new CanonicalBranchGradientCache(
                timeCount, stateDimension, transitionRepresentation, timeGrid);
        this.smootherStats = new BranchSmootherStats[timeCount];
        this.messageWorkspace = new CanonicalGaussianMessageOps.Workspace(stateDimension);
        for (int t = 0; t < timeCount; ++t) {
            smootherStats[t] = new BranchSmootherStats(t, stateDimension, t < timeCount - 1);
        }

        this.designMatrix = new double[observationDimension * stateDimension];
        this.noiseCovariance = new double[observationDimension * observationDimension];
        this.noisePrecision = new double[observationDimension * observationDimension];
        this.observationPrecisionContribution = new double[stateDimension * stateDimension];
        this.obsWorkspace = new double[observationDimension * stateDimension];
        this.observationVector = new double[observationDimension];
        this.observationInformation = new double[stateDimension];
        this.observationVectorWorkspace = new double[observationDimension];
        this.stateVectorWorkspace = new double[stateDimension];
        this.flatMatrixWorkspace = new double[maximumMatrixDimension * maximumMatrixDimension];
        this.flatInverseWorkspace = new double[maximumMatrixDimension * maximumMatrixDimension];
        this.flatCholeskyWorkspace = new double[maximumMatrixDimension * maximumMatrixDimension];
        this.flatLowerInverseWorkspace = new double[maximumMatrixDimension * maximumMatrixDimension];
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

    public double[][] getSmoothedMeans() {
        ensureResults();
        ensureMomentTrajectory();
        final double[][] out = new double[timeCount][stateDimension];
        for (int t = 0; t < timeCount; ++t) {
            System.arraycopy(smootherStats[t].smoothedMean, 0, out[t], 0, stateDimension);
        }
        return out;
    }

    public double[][][] getSmoothedCovariances() {
        ensureResults();
        ensureMomentTrajectory();
        final double[][][] out = new double[timeCount][stateDimension][stateDimension];
        for (int t = 0; t < timeCount; ++t) {
            MatrixOps.fromFlat(smootherStats[t].smoothedCovariance, out[t], stateDimension);
        }
        return out;
    }

    public double[][] getFilteredMeans() {
        ensureResults();
        ensureMomentTrajectory();
        final double[][] out = new double[timeCount][stateDimension];
        for (int t = 0; t < timeCount; ++t) {
            trajectory.copyFilteredMeanTo(t, out[t]);
        }
        return out;
    }

    public double[][][] getFilteredCovariances() {
        ensureResults();
        ensureMomentTrajectory();
        final double[][][] out = new double[timeCount][stateDimension][stateDimension];
        for (int t = 0; t < timeCount; ++t) {
            trajectory.copyFilteredCovarianceTo(t, out[t]);
        }
        return out;
    }

    public double[][] getPredictedMeans() {
        ensureResults();
        ensureMomentTrajectory();
        final double[][] out = new double[timeCount][stateDimension];
        for (int t = 0; t < timeCount; ++t) {
            trajectory.copyPredictedMeanTo(t, out[t]);
        }
        return out;
    }

    public double[][][] getPredictedCovariances() {
        ensureResults();
        ensureMomentTrajectory();
        final double[][][] out = new double[timeCount][stateDimension][stateDimension];
        for (int t = 0; t < timeCount; ++t) {
            trajectory.copyPredictedCovarianceTo(t, out[t]);
        }
        return out;
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
    public BranchSmootherStats[] getSmootherStats() {
        ensureResults();
        ensureMomentTrajectory();
        return smootherStats;
    }

    @Override
    public ForwardTrajectory getTrajectory() {
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
        observationModel.fillDesignMatrixFlat(designMatrix, stateDimension);
        observationModel.fillNoiseCovarianceFlat(noiseCovariance);
        final double noiseLogDet = invertPositiveDefinite(noiseCovariance, noisePrecision, observationDimension);
        buildObservationPrecisionContribution();

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

            if (observationModel.isObservationMissing(timeIndex)) {
                copyState(predictedState, filteredState);
            } else {
                observationModel.fillObservationVector(timeIndex, observationVector);
                buildObservationInformation(observationVector);

                addMatrices(predictedState.precision, observationPrecisionContribution, filteredState.precision);
                addVectors(predictedState.information, observationInformation, filteredState.information);
                filteredState.logNormalizer = normalizedLogNormalizerFlat(filteredState.precision,
                        filteredState.information);

                value += filteredState.logNormalizer
                        - predictedState.logNormalizer
                        - observationPotentialLogNormalizer(observationVector, noiseLogDet);
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

    private void buildObservationPrecisionContribution() {
        MatrixOps.matMul(noisePrecision, designMatrix, obsWorkspace,
                observationDimension, observationDimension, stateDimension);
        MatrixOps.matMulTransposedLeft(designMatrix, obsWorkspace, observationPrecisionContribution,
                observationDimension, stateDimension);
        MatrixOps.symmetrize(observationPrecisionContribution, stateDimension);
    }

    private void buildObservationInformation(final double[] observation) {
        MatrixOps.matVec(noisePrecision, observation, observationVectorWorkspace,
                observationDimension, observationDimension);
        for (int i = 0; i < stateDimension; ++i) {
            double sum = 0.0;
            for (int j = 0; j < observationDimension; ++j) {
                sum += designMatrix[j * stateDimension + i] * observationVectorWorkspace[j];
            }
            observationInformation[i] = sum;
        }
    }

    private double observationPotentialLogNormalizer(final double[] observation, final double logDetNoise) {
        final double quadratic = MatrixOps.quadraticForm(
                observation, noisePrecision, observationDimension, observationVectorWorkspace);
        return 0.5 * (observationDimension * MatrixOps.LOG_TWO_PI + logDetNoise + quadratic);
    }

    private double normalizedLogNormalizerFlat(final double[] precision, final double[] information) {
        if (!MatrixOps.tryCholesky(precision, flatCholeskyWorkspace, stateDimension)) {
            throw new IllegalArgumentException("Matrix is not positive definite");
        }
        final double logDet = MatrixOps.invertFromCholesky(
                flatCholeskyWorkspace, flatLowerInverseWorkspace, flatInverseWorkspace, stateDimension);
        MatrixOps.matVec(flatInverseWorkspace, information, stateVectorWorkspace, stateDimension);
        final double quadratic = dot(information, stateVectorWorkspace);
        return 0.5 * (stateDimension * MatrixOps.LOG_TWO_PI - logDet + quadratic);
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

    private double invertPositiveDefinite(final double[] matrix,
                                          final double[] inverseOut,
                                          final int dimension) {
        if (!MatrixOps.tryCholesky(matrix, flatCholeskyWorkspace, dimension)) {
            throw new IllegalArgumentException("Matrix is not positive definite");
        }
        return MatrixOps.invertFromCholesky(
                flatCholeskyWorkspace, flatLowerInverseWorkspace, inverseOut, dimension);
    }

    private double validatedDelta(final int fromIndex, final int toIndex) {
        final double dt = timeGrid.getDelta(fromIndex, toIndex);
        if (!(dt > 0.0)) {
            throw new IllegalArgumentException("Time increments must be strictly positive");
        }
        return dt;
    }

    private static void addMatrices(final double[] left, final double[] right, final double[] out) {
        for (int i = 0; i < left.length; ++i) {
            out[i] = left[i] + right[i];
        }
    }

    private static void addVectors(final double[] left, final double[] right, final double[] out) {
        for (int i = 0; i < left.length; ++i) {
            out[i] = left[i] + right[i];
        }
    }

    private static double dot(final double[] left, final double[] right) {
        double sum = 0.0;
        for (int i = 0; i < left.length; ++i) {
            sum += left[i] * right[i];
        }
        return sum;
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
