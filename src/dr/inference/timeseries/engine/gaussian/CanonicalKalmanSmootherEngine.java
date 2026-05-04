package dr.inference.timeseries.engine.gaussian;

import dr.evomodel.treedatalikelihood.continuous.canonical.message.GaussianMatrixOps;

import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianMessageOps;

import dr.evomodel.treedatalikelihood.continuous.canonical.math.MatrixOps;
import dr.inference.timeseries.core.TimeGrid;
import dr.inference.timeseries.gaussian.GaussianObservationModel;
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
    private final GaussianObservationModel observationModel;
    private final TimeGrid timeGrid;

    private final int stateDimension;
    private final int timeCount;
    private final int observationDimension;
    private final int maximumMatrixDimension;

    private final CanonicalGaussianState filteredCanonical;
    private final CanonicalGaussianState predictedCanonical;
    private final CanonicalGaussianTransition transitionCanonical;
    private final CanonicalGaussianState futureMessage;
    private final CanonicalGaussianState backwardMessage;
    private final ForwardTrajectory trajectory;
    private final CanonicalForwardTrajectory canonicalTrajectory;
    private final CanonicalBranchGradientCache branchGradientCache;
    private final BranchSmootherStats[] smootherStats;
    private final CanonicalGaussianMessageOps.Workspace messageWorkspace;

    private final double[][] designMatrix;
    private final double[][] noiseCovariance;
    private final double[][] noisePrecision;
    private final double[][] observationPrecisionContribution;
    private final double[][] transitionWorkspace;
    private final double[][] stateWorkspace;
    private final double[][] stateWorkspace2;
    private final double[][] stateWorkspace3;
    private final double[][] obsWorkspace;
    private final double[][] tempDxD1;
    private final double[][] tempDxD2;
    private final double[] observationVector;
    private final double[] observationInformation;
    private final double[] observationVectorWorkspace;
    private final double[] stateVectorWorkspace;
    private final double[] stateVectorWorkspace2;
    private final double[] tempD;
    private final double[] flatMatrixWorkspace;
    private final double[] flatInverseWorkspace;
    private final double[] flatCholeskyWorkspace;
    private final double[] flatLowerInverseWorkspace;

    private boolean resultsKnown;
    private boolean momentTrajectoryKnown;
    private double logLikelihood;

    public CanonicalKalmanSmootherEngine(final CanonicalGaussianBranchTransitionKernel canonicalKernel,
                                         final GaussianTransitionRepresentation transitionRepresentation,
                                         final GaussianObservationModel observationModel,
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

        this.filteredCanonical = new CanonicalGaussianState(stateDimension);
        this.predictedCanonical = new CanonicalGaussianState(stateDimension);
        this.transitionCanonical = new CanonicalGaussianTransition(stateDimension);

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

        this.designMatrix = new double[observationDimension][stateDimension];
        this.noiseCovariance = new double[observationDimension][observationDimension];
        this.noisePrecision = new double[observationDimension][observationDimension];
        this.observationPrecisionContribution = new double[stateDimension][stateDimension];
        this.transitionWorkspace = new double[stateDimension][stateDimension];
        this.stateWorkspace = new double[stateDimension][stateDimension];
        this.stateWorkspace2 = new double[stateDimension][stateDimension];
        this.stateWorkspace3 = new double[stateDimension][stateDimension];
        this.obsWorkspace = new double[observationDimension][stateDimension];
        this.tempDxD1 = new double[stateDimension][stateDimension];
        this.tempDxD2 = new double[stateDimension][stateDimension];
        this.observationVector = new double[observationDimension];
        this.observationInformation = new double[stateDimension];
        this.observationVectorWorkspace = new double[observationDimension];
        this.stateVectorWorkspace = new double[stateDimension];
        this.stateVectorWorkspace2 = new double[stateDimension];
        this.tempD = new double[stateDimension];
        this.flatMatrixWorkspace = new double[maximumMatrixDimension * maximumMatrixDimension];
        this.flatInverseWorkspace = new double[maximumMatrixDimension * maximumMatrixDimension];
        this.flatCholeskyWorkspace = new double[maximumMatrixDimension * maximumMatrixDimension];
        this.flatLowerInverseWorkspace = new double[maximumMatrixDimension * maximumMatrixDimension];
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
    }

    public double[][] getSmoothedMeans() {
        ensureResults();
        ensureMomentTrajectory();
        final double[][] out = new double[timeCount][stateDimension];
        for (int t = 0; t < timeCount; ++t) {
            GaussianMatrixOps.copyVector(smootherStats[t].smoothedMean, out[t]);
        }
        return out;
    }

    public double[][][] getSmoothedCovariances() {
        ensureResults();
        ensureMomentTrajectory();
        final double[][][] out = new double[timeCount][stateDimension][stateDimension];
        for (int t = 0; t < timeCount; ++t) {
            GaussianMatrixOps.copyMatrix(smootherStats[t].smoothedCovariance, out[t]);
        }
        return out;
    }

    public double[][] getFilteredMeans() {
        ensureResults();
        ensureMomentTrajectory();
        final double[][] out = new double[timeCount][stateDimension];
        for (int t = 0; t < timeCount; ++t) {
            GaussianMatrixOps.copyVector(trajectory.filteredMeans[t], out[t]);
        }
        return out;
    }

    public double[][][] getFilteredCovariances() {
        ensureResults();
        ensureMomentTrajectory();
        final double[][][] out = new double[timeCount][stateDimension][stateDimension];
        for (int t = 0; t < timeCount; ++t) {
            GaussianMatrixOps.copyMatrix(trajectory.filteredCovariances[t], out[t]);
        }
        return out;
    }

    public double[][] getPredictedMeans() {
        ensureResults();
        ensureMomentTrajectory();
        final double[][] out = new double[timeCount][stateDimension];
        for (int t = 0; t < timeCount; ++t) {
            GaussianMatrixOps.copyVector(trajectory.predictedMeans[t], out[t]);
        }
        return out;
    }

    public double[][][] getPredictedCovariances() {
        ensureResults();
        ensureMomentTrajectory();
        final double[][][] out = new double[timeCount][stateDimension][stateDimension];
        for (int t = 0; t < timeCount; ++t) {
            GaussianMatrixOps.copyMatrix(trajectory.predictedCovariances[t], out[t]);
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

    private double runForwardPass() {
        observationModel.fillDesignMatrix(designMatrix);
        observationModel.fillNoiseCovariance(noiseCovariance);
        final double noiseLogDet = invertPositiveDefinite(noiseCovariance, noisePrecision, observationDimension);
        buildObservationPrecisionContribution();

        canonicalKernel.fillInitialCanonicalState(filteredCanonical);

        double value = 0.0;
        for (int timeIndex = 0; timeIndex < timeCount; ++timeIndex) {
            if (timeIndex == 0) {
                copyState(filteredCanonical, predictedCanonical);
            } else {
                fillCanonicalTransition(timeIndex - 1, timeIndex, transitionCanonical);
                predict(filteredCanonical, transitionCanonical, predictedCanonical);
                copyTransition(transitionCanonical, canonicalTrajectory.transitions[timeIndex - 1]);

            }

            copyState(predictedCanonical, canonicalTrajectory.predictedStates[timeIndex]);

            if (observationModel.isObservationMissing(timeIndex)) {
                copyState(predictedCanonical, filteredCanonical);
            } else {
                observationModel.fillObservationVector(timeIndex, observationVector);
                buildObservationInformation(observationVector);

                addMatricesFlatAndRagged(predictedCanonical.precision, observationPrecisionContribution, filteredCanonical.precision, stateDimension);
                addVectors(predictedCanonical.information, observationInformation, filteredCanonical.information);
                filteredCanonical.logNormalizer = normalizedLogNormalizerFlat(filteredCanonical.precision,
                        filteredCanonical.information);

                value += filteredCanonical.logNormalizer
                        - predictedCanonical.logNormalizer
                        - observationPotentialLogNormalizer(observationVector, noiseLogDet);
            }

            copyState(filteredCanonical, canonicalTrajectory.filteredStates[timeIndex]);
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
                    trajectory.predictedMeans[timeIndex],
                    trajectory.predictedCovariances[timeIndex]);
            fillMomentsFromCanonical(canonicalTrajectory.filteredStates[timeIndex],
                    trajectory.filteredMeans[timeIndex],
                    trajectory.filteredCovariances[timeIndex]);
            fillMomentsFromCanonical(canonicalTrajectory.smoothedStates[timeIndex],
                    smootherStats[timeIndex].smoothedMean,
                    smootherStats[timeIndex].smoothedCovariance);
            if (timeIndex < timeCount - 1) {
                transitionRepresentation.getTransitionMatrix(timeIndex, timeIndex + 1, timeGrid,
                        trajectory.transitionMatrices[timeIndex]);
                transitionRepresentation.getTransitionOffset(timeIndex, timeIndex + 1, timeGrid,
                        trajectory.transitionOffsets[timeIndex]);
                transitionRepresentation.getTransitionCovariance(timeIndex, timeIndex + 1, timeGrid,
                        trajectory.stepCovariances[timeIndex]);
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
                                          final double[][] covarianceOut) {
        invertPositiveDefiniteFlatInput(canonical.precision, covarianceOut, stateDimension);
        GaussianMatrixOps.multiplyMatrixVector(covarianceOut, canonical.information, meanOut,
                stateDimension, stateDimension);
    }

    private void buildObservationPrecisionContribution() {
        GaussianMatrixOps.multiplyMatrixMatrix(noisePrecision, designMatrix, obsWorkspace,
                observationDimension, observationDimension, stateDimension);
        GaussianMatrixOps.multiplyMatrixMatrixTransposedRight(
                designMatrix, obsWorkspace, observationPrecisionContribution);
        GaussianMatrixOps.symmetrize(observationPrecisionContribution);
    }

    private void buildObservationInformation(final double[] observation) {
        GaussianMatrixOps.multiplyMatrixVector(noisePrecision, observation, observationVectorWorkspace,
                observationDimension, observationDimension);
        for (int i = 0; i < stateDimension; ++i) {
            double sum = 0.0;
            for (int j = 0; j < observationDimension; ++j) {
                sum += designMatrix[j][i] * observationVectorWorkspace[j];
            }
            observationInformation[i] = sum;
        }
    }

    private double observationPotentialLogNormalizer(final double[] observation, final double logDetNoise) {
        final double quadratic = GaussianMatrixOps.quadraticForm(noisePrecision, observation);
        return 0.5 * (observationDimension * GaussianMatrixOps.LOG_TWO_PI + logDetNoise + quadratic);
    }

    private double normalizedLogNormalizer(final double[][] precision, final double[] information) {
        final double[][] precisionInverse = stateWorkspace;
        final double logDet = invertPositiveDefinite(precision, precisionInverse, stateDimension);
        GaussianMatrixOps.multiplyMatrixVector(precisionInverse, information, stateVectorWorkspace,
                stateDimension, stateDimension);
        final double quadratic = dot(information, stateVectorWorkspace);
        return 0.5 * (stateDimension * GaussianMatrixOps.LOG_TWO_PI - logDet + quadratic);
    }

    private double normalizedLogNormalizerFlat(final double[] precision, final double[] information) {
        if (!MatrixOps.tryCholesky(precision, flatCholeskyWorkspace, stateDimension)) {
            throw new IllegalArgumentException("Matrix is not positive definite");
        }
        final double logDet = MatrixOps.invertFromCholesky(
                flatCholeskyWorkspace, flatLowerInverseWorkspace, flatInverseWorkspace, stateDimension);
        MatrixOps.matVec(flatInverseWorkspace, information, stateVectorWorkspace, stateDimension);
        final double quadratic = dot(information, stateVectorWorkspace);
        return 0.5 * (stateDimension * GaussianMatrixOps.LOG_TWO_PI - logDet + quadratic);
    }

    private void invertPositiveDefiniteFlatInput(final double[] matrix,
                                                  final double[][] inverseOut,
                                                  final int dimension) {
        if (!MatrixOps.tryCholesky(matrix, flatCholeskyWorkspace, dimension)) {
            throw new IllegalArgumentException("Matrix is not positive definite");
        }
        MatrixOps.invertFromCholesky(flatCholeskyWorkspace, flatLowerInverseWorkspace, flatInverseWorkspace, dimension);
        copyFlatToMatrix(flatInverseWorkspace, inverseOut, dimension);
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

    private static void addMatricesFlatAndRagged(final double[] left, final double[][] right,
                                                  final double[] out, final int dim) {
        for (int i = 0; i < dim; ++i) {
            final int rowOff = i * dim;
            for (int j = 0; j < dim; ++j) {
                out[rowOff + j] = left[rowOff + j] + right[i][j];
            }
        }
    }

    private double invertPositiveDefinite(final double[][] matrix,
                                          final double[][] inverseOut,
                                          final int dimension) {
        copyMatrixToFlat(matrix, flatMatrixWorkspace, dimension);
        if (!MatrixOps.tryCholesky(flatMatrixWorkspace, flatCholeskyWorkspace, dimension)) {
            throw new IllegalArgumentException("Matrix is not positive definite");
        }
        final double logDet = MatrixOps.invertFromCholesky(
                flatCholeskyWorkspace, flatLowerInverseWorkspace, flatInverseWorkspace, dimension);
        copyFlatToMatrix(flatInverseWorkspace, inverseOut, dimension);
        return logDet;
    }

    private static void copyMatrixToFlat(final double[][] source,
                                         final double[] target,
                                         final int dimension) {
        for (int i = 0; i < dimension; ++i) {
            System.arraycopy(source[i], 0, target, i * dimension, dimension);
        }
    }

    private static void copyFlatToMatrix(final double[] source,
                                         final double[][] target,
                                         final int dimension) {
        for (int i = 0; i < dimension; ++i) {
            System.arraycopy(source, i * dimension, target[i], 0, dimension);
        }
    }

    private double validatedDelta(final int fromIndex, final int toIndex) {
        final double dt = timeGrid.getDelta(fromIndex, toIndex);
        if (!(dt > 0.0)) {
            throw new IllegalArgumentException("Time increments must be strictly positive");
        }
        return dt;
    }

    private static void addMatrices(final double[][] left, final double[][] right, final double[][] out) {
        for (int i = 0; i < left.length; ++i) {
            for (int j = 0; j < left[i].length; ++j) {
                out[i][j] = left[i][j] + right[i][j];
            }
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
        GaussianMatrixOps.copyVector(source.information, target.information);
        target.logNormalizer = source.logNormalizer;
    }

    private static void copyTransition(final CanonicalGaussianTransition source,
                                       final CanonicalGaussianTransition target) {
        System.arraycopy(source.precisionXX, 0, target.precisionXX, 0, source.precisionXX.length);
        System.arraycopy(source.precisionXY, 0, target.precisionXY, 0, source.precisionXY.length);
        System.arraycopy(source.precisionYX, 0, target.precisionYX, 0, source.precisionYX.length);
        System.arraycopy(source.precisionYY, 0, target.precisionYY, 0, source.precisionYY.length);
        GaussianMatrixOps.copyVector(source.informationX, target.informationX);
        GaussianMatrixOps.copyVector(source.informationY, target.informationY);
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

    private void fillUpperLeftBlock(final double[][] source, final double[][] out) {
        for (int i = 0; i < stateDimension; ++i) {
            System.arraycopy(source[i], 0, out[i], 0, stateDimension);
        }
    }

    private void fillUpperRightBlock(final double[][] source, final double[][] out) {
        for (int i = 0; i < stateDimension; ++i) {
            for (int j = 0; j < stateDimension; ++j) {
                out[i][j] = source[i][stateDimension + j];
            }
        }
    }

    private void fillLowerLeftBlock(final double[][] source, final double[][] out) {
        for (int i = 0; i < stateDimension; ++i) {
            for (int j = 0; j < stateDimension; ++j) {
                out[i][j] = source[stateDimension + i][j];
            }
        }
    }

    private void fillLowerRightBlock(final double[][] source, final double[][] out) {
        for (int i = 0; i < stateDimension; ++i) {
            for (int j = 0; j < stateDimension; ++j) {
                out[i][j] = source[stateDimension + i][stateDimension + j];
            }
        }
    }

    private void fillFirstBlock(final double[] source, final double[] out) {
        System.arraycopy(source, 0, out, 0, stateDimension);
    }

    private void fillSecondBlock(final double[] source, final double[] out) {
        System.arraycopy(source, stateDimension, out, 0, stateDimension);
    }

    private void multiply(final double[][] left, final double[][] right, final double[][] out) {
        for (int i = 0; i < left.length; ++i) {
            for (int j = 0; j < right[0].length; ++j) {
                double sum = 0.0;
                for (int k = 0; k < right.length; ++k) {
                    sum += left[i][k] * right[k][j];
                }
                out[i][j] = sum;
            }
        }
    }

    private void subtractMatrixInPlace(final double[][] target, final double[][] delta) {
        for (int i = 0; i < target.length; ++i) {
            for (int j = 0; j < target[i].length; ++j) {
                target[i][j] -= delta[i][j];
            }
        }
    }

    private double normalizedLogNormalizer(final double[][] precision,
                                           final double[] information,
                                           final int dimension) {
        final double[][] precisionInverse = stateWorkspace;
        final double logDet = invertPositiveDefinite(precision, precisionInverse, dimension);
        GaussianMatrixOps.multiplyMatrixVector(precisionInverse, information, stateVectorWorkspace,
                dimension, dimension);
        final double quadratic = dot(information, stateVectorWorkspace);
        return 0.5 * (dimension * GaussianMatrixOps.LOG_TWO_PI - logDet + quadratic);
    }
}
