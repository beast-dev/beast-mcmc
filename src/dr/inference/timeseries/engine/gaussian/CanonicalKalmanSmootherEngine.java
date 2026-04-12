package dr.inference.timeseries.engine.gaussian;

import dr.inference.timeseries.core.TimeGrid;
import dr.inference.timeseries.gaussian.GaussianObservationModel;
import dr.inference.timeseries.representation.CanonicalGaussianBranchTransitionKernel;
import dr.inference.timeseries.representation.CanonicalGaussianState;
import dr.inference.timeseries.representation.CanonicalGaussianTransition;
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
    private final GaussianObservationModel observationModel;
    private final TimeGrid timeGrid;

    private final int stateDimension;
    private final int timeCount;
    private final int observationDimension;

    private final CanonicalGaussianState filteredCanonical;
    private final CanonicalGaussianState predictedCanonical;
    private final CanonicalGaussianTransition transitionCanonical;
    private final CanonicalGaussianState futureMessage;
    private final ForwardTrajectory trajectory;
    private final CanonicalForwardTrajectory canonicalTrajectory;
    private final BranchSmootherStats[] smootherStats;

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

    private boolean resultsKnown;
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
        this.observationModel = observationModel;
        this.timeGrid = timeGrid;
        this.stateDimension = canonicalKernel.getStateDimension();
        this.timeCount = timeGrid.getTimeCount();
        this.observationDimension = observationModel.getObservationDimension();

        this.filteredCanonical = new CanonicalGaussianState(stateDimension);
        this.predictedCanonical = new CanonicalGaussianState(stateDimension);
        this.transitionCanonical = new CanonicalGaussianTransition(stateDimension);

        this.futureMessage = new CanonicalGaussianState(stateDimension);
        this.trajectory = new ForwardTrajectory(timeCount, stateDimension);
        this.canonicalTrajectory = new CanonicalForwardTrajectory(timeCount, stateDimension);
        this.smootherStats = new BranchSmootherStats[timeCount];
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
    }

    public double getLogLikelihood() {
        ensureResults();
        return logLikelihood;
    }

    @Override
    public void makeDirty() {
        resultsKnown = false;
    }

    public double[][] getSmoothedMeans() {
        ensureResults();
        final double[][] out = new double[timeCount][stateDimension];
        for (int t = 0; t < timeCount; ++t) {
            KalmanLikelihoodEngine.copyVector(smootherStats[t].smoothedMean, out[t]);
        }
        return out;
    }

    public double[][][] getSmoothedCovariances() {
        ensureResults();
        final double[][][] out = new double[timeCount][stateDimension][stateDimension];
        for (int t = 0; t < timeCount; ++t) {
            KalmanLikelihoodEngine.copyMatrix(smootherStats[t].smoothedCovariance, out[t]);
        }
        return out;
    }

    public double[][] getFilteredMeans() {
        ensureResults();
        final double[][] out = new double[timeCount][stateDimension];
        for (int t = 0; t < timeCount; ++t) {
            KalmanLikelihoodEngine.copyVector(trajectory.filteredMeans[t], out[t]);
        }
        return out;
    }

    public double[][][] getFilteredCovariances() {
        ensureResults();
        final double[][][] out = new double[timeCount][stateDimension][stateDimension];
        for (int t = 0; t < timeCount; ++t) {
            KalmanLikelihoodEngine.copyMatrix(trajectory.filteredCovariances[t], out[t]);
        }
        return out;
    }

    public double[][] getPredictedMeans() {
        ensureResults();
        final double[][] out = new double[timeCount][stateDimension];
        for (int t = 0; t < timeCount; ++t) {
            KalmanLikelihoodEngine.copyVector(trajectory.predictedMeans[t], out[t]);
        }
        return out;
    }

    public double[][][] getPredictedCovariances() {
        ensureResults();
        final double[][][] out = new double[timeCount][stateDimension][stateDimension];
        for (int t = 0; t < timeCount; ++t) {
            KalmanLikelihoodEngine.copyMatrix(trajectory.predictedCovariances[t], out[t]);
        }
        return out;
    }

    private void ensureResults() {
        if (!resultsKnown) {
            logLikelihood = runForwardPass();
            runBackwardPass();
            resultsKnown = true;
        }
    }

    @Override
    public BranchSmootherStats[] getSmootherStats() {
        ensureResults();
        return smootherStats;
    }

    @Override
    public ForwardTrajectory getTrajectory() {
        ensureResults();
        return trajectory;
    }

    @Override
    public CanonicalForwardTrajectory getCanonicalTrajectory() {
        ensureResults();
        return canonicalTrajectory;
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
                final double dt = validatedDelta(timeIndex - 1, timeIndex);
                canonicalKernel.fillCanonicalTransition(dt, transitionCanonical);
                predict(filteredCanonical, transitionCanonical, predictedCanonical);
                copyTransition(transitionCanonical, canonicalTrajectory.transitions[timeIndex - 1]);

                transitionRepresentation.getTransitionMatrix(timeIndex - 1, timeIndex, timeGrid,
                        trajectory.transitionMatrices[timeIndex - 1]);
                transitionRepresentation.getTransitionOffset(timeIndex - 1, timeIndex, timeGrid,
                        trajectory.transitionOffsets[timeIndex - 1]);
                transitionRepresentation.getTransitionCovariance(timeIndex - 1, timeIndex, timeGrid,
                        trajectory.stepCovariances[timeIndex - 1]);
            }

            fillMomentsFromCanonical(predictedCanonical,
                    trajectory.predictedMeans[timeIndex],
                    trajectory.predictedCovariances[timeIndex]);
            copyState(predictedCanonical, canonicalTrajectory.predictedStates[timeIndex]);

            if (observationModel.isObservationMissing(timeIndex)) {
                copyState(predictedCanonical, filteredCanonical);
            } else {
                observationModel.fillObservationVector(timeIndex, observationVector);
                buildObservationInformation(observationVector);

                addMatrices(predictedCanonical.precision, observationPrecisionContribution, filteredCanonical.precision);
                addVectors(predictedCanonical.information, observationInformation, filteredCanonical.information);
                filteredCanonical.logNormalizer = normalizedLogNormalizer(filteredCanonical.precision,
                        filteredCanonical.information);

                value += filteredCanonical.logNormalizer
                        - predictedCanonical.logNormalizer
                        - observationPotentialLogNormalizer(observationVector, noiseLogDet);
            }

            fillMomentsFromCanonical(filteredCanonical,
                    trajectory.filteredMeans[timeIndex],
                    trajectory.filteredCovariances[timeIndex]);
            copyState(filteredCanonical, canonicalTrajectory.filteredStates[timeIndex]);
        }
        return value;
    }

    private void runBackwardPass() {
        final int T = trajectory.timeCount;
        copyState(canonicalTrajectory.filteredStates[T - 1], canonicalTrajectory.smoothedStates[T - 1]);
        fillMomentsFromCanonical(canonicalTrajectory.smoothedStates[T - 1],
                smootherStats[T - 1].smoothedMean,
                smootherStats[T - 1].smoothedCovariance);

        for (int t = T - 2; t >= 0; --t) {
            subtractStates(canonicalTrajectory.smoothedStates[t + 1],
                    canonicalTrajectory.predictedStates[t + 1],
                    futureMessage);
            buildPairPosterior(canonicalTrajectory.filteredStates[t],
                    canonicalTrajectory.transitions[t],
                    futureMessage,
                    canonicalTrajectory.branchPairStates[t]);
            marginalizeFirstBlock(canonicalTrajectory.branchPairStates[t], canonicalTrajectory.smoothedStates[t]);
            fillMomentsFromCanonical(canonicalTrajectory.smoothedStates[t],
                    smootherStats[t].smoothedMean,
                    smootherStats[t].smoothedCovariance);
        }
    }

    private void predict(final CanonicalGaussianState previous,
                         final CanonicalGaussianTransition transition,
                         final CanonicalGaussianState out) {
        addMatrices(previous.precision, transition.precisionXX, stateWorkspace);
        addVectors(previous.information, transition.informationX, stateVectorWorkspace);

        final double[][] statePrecisionInverse = stateWorkspace2;
        invertPositiveDefinite(stateWorkspace, statePrecisionInverse, stateDimension);

        KalmanLikelihoodEngine.multiplyMatrixMatrix(statePrecisionInverse, transition.precisionXY, transitionWorkspace);
        KalmanLikelihoodEngine.multiplyMatrixMatrix(transition.precisionYX, transitionWorkspace, stateWorkspace3);

        for (int i = 0; i < stateDimension; ++i) {
            for (int j = 0; j < stateDimension; ++j) {
                out.precision[i][j] = transition.precisionYY[i][j] - stateWorkspace3[i][j];
            }
        }
        KalmanLikelihoodEngine.symmetrize(out.precision);

        KalmanLikelihoodEngine.multiplyMatrixVector(statePrecisionInverse, stateVectorWorkspace, stateVectorWorkspace2,
                stateDimension, stateDimension);
        KalmanLikelihoodEngine.multiplyMatrixVector(transition.precisionYX, stateVectorWorkspace2, stateVectorWorkspace,
                stateDimension, stateDimension);
        for (int i = 0; i < stateDimension; ++i) {
            out.information[i] = transition.informationY[i] - stateVectorWorkspace[i];
        }

        out.logNormalizer = normalizedLogNormalizer(out.precision, out.information);
    }

    private void fillMomentsFromCanonical(final CanonicalGaussianState canonical,
                                          final double[] meanOut,
                                          final double[][] covarianceOut) {
        invertPositiveDefinite(canonical.precision, covarianceOut, stateDimension);
        KalmanLikelihoodEngine.multiplyMatrixVector(covarianceOut, canonical.information, meanOut,
                stateDimension, stateDimension);
    }

    private void buildObservationPrecisionContribution() {
        KalmanLikelihoodEngine.multiplyMatrixMatrix(noisePrecision, designMatrix, obsWorkspace,
                observationDimension, observationDimension, stateDimension);
        KalmanLikelihoodEngine.multiplyMatrixMatrixTransposedRight(
                designMatrix, obsWorkspace, observationPrecisionContribution);
        KalmanLikelihoodEngine.symmetrize(observationPrecisionContribution);
    }

    private void buildObservationInformation(final double[] observation) {
        KalmanLikelihoodEngine.multiplyMatrixVector(noisePrecision, observation, observationVectorWorkspace,
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
        final double quadratic = KalmanLikelihoodEngine.quadraticForm(noisePrecision, observation);
        return 0.5 * (observationDimension * KalmanLikelihoodEngine.LOG_TWO_PI + logDetNoise + quadratic);
    }

    private double normalizedLogNormalizer(final double[][] precision, final double[] information) {
        final double[][] precisionInverse = stateWorkspace;
        final double logDet = invertPositiveDefinite(precision, precisionInverse, stateDimension);
        KalmanLikelihoodEngine.multiplyMatrixVector(precisionInverse, information, stateVectorWorkspace,
                stateDimension, stateDimension);
        final double quadratic = dot(information, stateVectorWorkspace);
        return 0.5 * (stateDimension * KalmanLikelihoodEngine.LOG_TWO_PI - logDet + quadratic);
    }

    private static double invertPositiveDefinite(final double[][] matrix,
                                                 final double[][] inverseOut,
                                                 final int dimension) {
        final double[][] copy = new double[dimension][dimension];
        KalmanLikelihoodEngine.copyMatrix(matrix, copy, dimension, dimension);
        final KalmanLikelihoodEngine.CholeskyFactor chol = KalmanLikelihoodEngine.cholesky(copy);
        KalmanLikelihoodEngine.copyMatrix(copy, inverseOut, dimension, dimension);
        KalmanLikelihoodEngine.invertPositiveDefiniteFromCholesky(inverseOut, chol);
        return chol.logDeterminant();
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
        KalmanLikelihoodEngine.copyMatrix(source.precision, target.precision);
        KalmanLikelihoodEngine.copyVector(source.information, target.information);
        target.logNormalizer = source.logNormalizer;
    }

    private static void copyTransition(final CanonicalGaussianTransition source,
                                       final CanonicalGaussianTransition target) {
        KalmanLikelihoodEngine.copyMatrix(source.precisionXX, target.precisionXX);
        KalmanLikelihoodEngine.copyMatrix(source.precisionXY, target.precisionXY);
        KalmanLikelihoodEngine.copyMatrix(source.precisionYX, target.precisionYX);
        KalmanLikelihoodEngine.copyMatrix(source.precisionYY, target.precisionYY);
        KalmanLikelihoodEngine.copyVector(source.informationX, target.informationX);
        KalmanLikelihoodEngine.copyVector(source.informationY, target.informationY);
        target.logNormalizer = source.logNormalizer;
    }

    private void subtractStates(final CanonicalGaussianState left,
                                final CanonicalGaussianState right,
                                final CanonicalGaussianState out) {
        for (int i = 0; i < stateDimension; ++i) {
            out.information[i] = left.information[i] - right.information[i];
            for (int j = 0; j < stateDimension; ++j) {
                out.precision[i][j] = left.precision[i][j] - right.precision[i][j];
            }
        }
        out.logNormalizer = left.logNormalizer - right.logNormalizer;
    }

    private void buildPairPosterior(final CanonicalGaussianState filteredState,
                                    final CanonicalGaussianTransition transition,
                                    final CanonicalGaussianState futureMessage,
                                    final CanonicalGaussianState pairOut) {
        final int d = stateDimension;
        for (int i = 0; i < d; ++i) {
            pairOut.information[i] = filteredState.information[i] + transition.informationX[i];
            pairOut.information[d + i] = futureMessage.information[i] + transition.informationY[i];
            for (int j = 0; j < d; ++j) {
                pairOut.precision[i][j] = filteredState.precision[i][j] + transition.precisionXX[i][j];
                pairOut.precision[i][d + j] = transition.precisionXY[i][j];
                pairOut.precision[d + i][j] = transition.precisionYX[i][j];
                pairOut.precision[d + i][d + j] = futureMessage.precision[i][j] + transition.precisionYY[i][j];
            }
        }
        pairOut.logNormalizer = 0.0;
    }

    private void marginalizeFirstBlock(final CanonicalGaussianState pairState,
                                       final CanonicalGaussianState out) {
        final int d = stateDimension;
        final double[][] yyInverse = tempDxD1;
        fillLowerRightBlock(pairState.precision, tempDxD2);
        invertPositiveDefinite(tempDxD2, yyInverse, d);

        fillUpperRightBlock(pairState.precision, transitionWorkspace);
        fillLowerLeftBlock(pairState.precision, stateWorkspace3);
        multiply(transitionWorkspace, yyInverse, stateWorkspace);
        KalmanLikelihoodEngine.multiplyMatrixMatrix(stateWorkspace, stateWorkspace3, stateWorkspace2);

        fillUpperLeftBlock(pairState.precision, out.precision);
        subtractMatrixInPlace(out.precision, stateWorkspace2);
        KalmanLikelihoodEngine.symmetrize(out.precision);

        fillSecondBlock(pairState.information, stateVectorWorkspace2);
        KalmanLikelihoodEngine.multiplyMatrixVector(yyInverse, stateVectorWorkspace2, stateVectorWorkspace,
                d, d);
        KalmanLikelihoodEngine.multiplyMatrixVector(transitionWorkspace, stateVectorWorkspace, tempD,
                d, d);
        fillFirstBlock(pairState.information, out.information);
        for (int i = 0; i < d; ++i) {
            out.information[i] -= tempD[i];
        }
        out.logNormalizer = normalizedLogNormalizer(out.precision, out.information);
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
        KalmanLikelihoodEngine.multiplyMatrixVector(precisionInverse, information, stateVectorWorkspace,
                dimension, dimension);
        final double quadratic = dot(information, stateVectorWorkspace);
        return 0.5 * (dimension * KalmanLikelihoodEngine.LOG_TWO_PI - logDet + quadratic);
    }
}
