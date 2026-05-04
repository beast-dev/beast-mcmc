package dr.inference.timeseries.engine.gaussian;

import dr.evomodel.treedatalikelihood.continuous.canonical.message.GaussianMatrixOps;

import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianMessageOps;

import dr.evomodel.treedatalikelihood.continuous.canonical.math.MatrixOps;
import dr.inference.timeseries.core.TimeGrid;
import dr.inference.timeseries.engine.LikelihoodEngine;
import dr.inference.timeseries.gaussian.GaussianObservationModel;
import dr.inference.timeseries.representation.CachedGaussianTransitionRepresentation;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianBranchTransitionKernel;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianState;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianTransition;

/**
 * Forward-only Gaussian likelihood engine that propagates canonical-form states.
 *
 * <p>This engine is intentionally parallel to {@link KalmanLikelihoodEngine}: it
 * computes the same marginal likelihood, but uses precision/information algebra
 * throughout the latent-state propagation. It is the first step toward a selectable
 * canonical backend without disturbing the working expectation-form smoother and
 * gradient code.
 */
public final class CanonicalKalmanLikelihoodEngine implements LikelihoodEngine {

    private final CanonicalGaussianBranchTransitionKernel transitionKernel;
    private final CachedGaussianTransitionRepresentation cachedTransitionRepresentation;
    private final GaussianObservationModel observationModel;
    private final TimeGrid timeGrid;

    private final int stateDimension;
    private final int observationDimension;
    private final int timeCount;
    private final int maximumMatrixDimension;

    private final CanonicalGaussianState filteredState;
    private final CanonicalGaussianState predictedState;
    private final CanonicalGaussianTransition transition;
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
    private final double[] observationVector;
    private final double[] observationInformation;
    private final double[] stateVectorWorkspace;
    private final double[] stateVectorWorkspace2;
    private final double[] observationVectorWorkspace;
    private final double[] flatMatrixWorkspace;
    private final double[] flatInverseWorkspace;
    private final double[] flatCholeskyWorkspace;
    private final double[] flatLowerInverseWorkspace;

    private boolean likelihoodKnown;
    private double logLikelihood;

    public CanonicalKalmanLikelihoodEngine(final CanonicalGaussianBranchTransitionKernel transitionKernel,
                                           final GaussianObservationModel observationModel,
                                           final TimeGrid timeGrid) {
        this(transitionKernel, null, observationModel, timeGrid);
    }

    public CanonicalKalmanLikelihoodEngine(final CanonicalGaussianBranchTransitionKernel transitionKernel,
                                           final CachedGaussianTransitionRepresentation cachedTransitionRepresentation,
                                           final GaussianObservationModel observationModel,
                                           final TimeGrid timeGrid) {
        if (transitionKernel == null) {
            throw new IllegalArgumentException("transitionKernel must not be null");
        }
        if (observationModel == null) {
            throw new IllegalArgumentException("observationModel must not be null");
        }
        if (timeGrid == null) {
            throw new IllegalArgumentException("timeGrid must not be null");
        }
        this.transitionKernel = transitionKernel;
        this.cachedTransitionRepresentation = cachedTransitionRepresentation;
        this.observationModel = observationModel;
        this.timeGrid = timeGrid;
        this.stateDimension = transitionKernel.getStateDimension();
        this.observationDimension = observationModel.getObservationDimension();
        this.timeCount = timeGrid.getTimeCount();
        this.maximumMatrixDimension = Math.max(stateDimension, observationDimension);

        if (observationModel.getTimeCount() != timeCount) {
            throw new IllegalArgumentException(
                    "Observation time count (" + observationModel.getTimeCount() +
                            ") must match time grid count (" + timeCount + ")");
        }

        this.filteredState = new CanonicalGaussianState(stateDimension);
        this.predictedState = new CanonicalGaussianState(stateDimension);
        this.transition = new CanonicalGaussianTransition(stateDimension);
        this.messageWorkspace = new CanonicalGaussianMessageOps.Workspace(stateDimension);

        this.designMatrix = new double[observationDimension][stateDimension];
        this.noiseCovariance = new double[observationDimension][observationDimension];
        this.noisePrecision = new double[observationDimension][observationDimension];
        this.observationPrecisionContribution = new double[stateDimension][stateDimension];
        this.transitionWorkspace = new double[stateDimension][stateDimension];
        this.stateWorkspace = new double[stateDimension][stateDimension];
        this.stateWorkspace2 = new double[stateDimension][stateDimension];
        this.stateWorkspace3 = new double[stateDimension][stateDimension];
        this.obsWorkspace = new double[observationDimension][stateDimension];
        this.observationVector = new double[observationDimension];
        this.observationInformation = new double[stateDimension];
        this.stateVectorWorkspace = new double[stateDimension];
        this.stateVectorWorkspace2 = new double[stateDimension];
        this.observationVectorWorkspace = new double[observationDimension];
        this.flatMatrixWorkspace = new double[maximumMatrixDimension * maximumMatrixDimension];
        this.flatInverseWorkspace = new double[maximumMatrixDimension * maximumMatrixDimension];
        this.flatCholeskyWorkspace = new double[maximumMatrixDimension * maximumMatrixDimension];
        this.flatLowerInverseWorkspace = new double[maximumMatrixDimension * maximumMatrixDimension];
        if (cachedTransitionRepresentation != null) {
            cachedTransitionRepresentation.prepareTimeGrid(timeGrid);
        }
    }

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

    private double computeLogLikelihood() {
        observationModel.fillDesignMatrix(designMatrix);
        observationModel.fillNoiseCovariance(noiseCovariance);
        final double noiseLogDet = invertPositiveDefinite(noiseCovariance, noisePrecision, observationDimension);
        buildObservationPrecisionContribution();

        transitionKernel.fillInitialCanonicalState(filteredState);

        double value = 0.0;
        for (int timeIndex = 0; timeIndex < timeCount; ++timeIndex) {
            if (timeIndex == 0) {
                copyState(filteredState, predictedState);
            } else {
                fillCanonicalTransition(timeIndex - 1, timeIndex, transition);
                predict(filteredState, transition, predictedState);
            }

            if (observationModel.isObservationMissing(timeIndex)) {
                copyState(predictedState, filteredState);
                continue;
            }

            observationModel.fillObservationVector(timeIndex, observationVector);
            buildObservationInformation(observationVector);

            addMatrices(predictedState.precision, observationPrecisionContribution, filteredState.precision);
            addVectors(predictedState.information, observationInformation, filteredState.information);
            filteredState.logNormalizer = CanonicalGaussianMessageOps.normalizedLogNormalizer(
                    filteredState, messageWorkspace);

            value += filteredState.logNormalizer
                    - predictedState.logNormalizer
                    - observationPotentialLogNormalizer(observationVector, noiseLogDet);
        }
        return value;
    }

    private void predict(final CanonicalGaussianState previous,
                         final CanonicalGaussianTransition transition,
                         final CanonicalGaussianState out) {
        CanonicalGaussianMessageOps.pushForward(previous, transition, messageWorkspace, out);
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

    private void fillCanonicalTransition(final int fromIndex,
                                         final int toIndex,
                                         final CanonicalGaussianTransition out) {
        if (cachedTransitionRepresentation != null) {
            cachedTransitionRepresentation.getCanonicalTransition(fromIndex, toIndex, timeGrid, out);
        } else {
            transitionKernel.fillCanonicalTransition(validatedDelta(fromIndex, toIndex), out);
        }
    }

    private static void addMatrices(final double[] left,
                                    final double[][] right,
                                    final double[] out) {
        final int rows = right.length;
        final int cols = right[0].length;
        for (int i = 0; i < rows; ++i) {
            final int rowOffset = i * cols;
            for (int j = 0; j < cols; ++j) {
                out[rowOffset + j] = left[rowOffset + j] + right[i][j];
            }
        }
    }

    private static void addVectors(final double[] left,
                                   final double[] right,
                                   final double[] out) {
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

    private static void copyState(final CanonicalGaussianState source,
                                  final CanonicalGaussianState target) {
        System.arraycopy(source.precision, 0, target.precision, 0, source.precision.length);
        GaussianMatrixOps.copyVector(source.information, target.information);
        target.logNormalizer = source.logNormalizer;
    }
}
