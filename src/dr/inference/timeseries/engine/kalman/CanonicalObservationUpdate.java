package dr.inference.timeseries.engine.kalman;

import dr.evomodel.treedatalikelihood.continuous.canonical.math.MatrixOps;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianMessageOps;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianState;
import dr.inference.timeseries.model.gaussian.LinearGaussianObservationModel;

/**
 * Canonical-form contribution from the observed coordinates at one time point.
 */
final class CanonicalObservationUpdate {

    private final LinearGaussianObservationModel observationModel;
    private final int stateDimension;
    private final int observationDimension;

    private final int[] observedIndices;
    private final double[] fullDesignMatrix;
    private final double[] fullNoiseCovariance;
    private final double[] observedDesignMatrix;
    private final double[] observedNoiseCovariance;
    private final double[] observedNoisePrecision;
    private final double[] observedVector;
    private final double[] noisePrecisionTimesDesign;
    private final double[] noisePrecisionTimesObservation;
    private final double[] observationPrecisionContribution;
    private final double[] observationInformationContribution;
    private final double[] observationWorkspace;
    private final double[] flatCholeskyWorkspace;
    private final double[] flatLowerInverseWorkspace;

    private int observedCount;
    private double potentialLogNormalizer;

    CanonicalObservationUpdate(final LinearGaussianObservationModel observationModel,
                               final int stateDimension) {
        this.observationModel = observationModel;
        this.stateDimension = stateDimension;
        this.observationDimension = observationModel.getObservationDimension();

        this.observedIndices = new int[observationDimension];
        this.fullDesignMatrix = new double[observationDimension * stateDimension];
        this.fullNoiseCovariance = new double[observationDimension * observationDimension];
        this.observedDesignMatrix = new double[observationDimension * stateDimension];
        this.observedNoiseCovariance = new double[observationDimension * observationDimension];
        this.observedNoisePrecision = new double[observationDimension * observationDimension];
        this.observedVector = new double[observationDimension];
        this.noisePrecisionTimesDesign = new double[observationDimension * stateDimension];
        this.noisePrecisionTimesObservation = new double[observationDimension];
        this.observationPrecisionContribution = new double[stateDimension * stateDimension];
        this.observationInformationContribution = new double[stateDimension];
        this.observationWorkspace = new double[observationDimension];
        this.flatCholeskyWorkspace = new double[observationDimension * observationDimension];
        this.flatLowerInverseWorkspace = new double[observationDimension * observationDimension];
    }

    void refreshStaticMatrices() {
        observationModel.fillDesignMatrixFlat(fullDesignMatrix, stateDimension);
        observationModel.fillNoiseCovarianceFlat(fullNoiseCovariance);
    }

    boolean prepareForTime(final int timeIndex) {
        if (timeIndex < 0 || timeIndex >= observationModel.getTimeCount()) {
            throw new IllegalArgumentException("timeIndex out of bounds: " + timeIndex);
        }

        observedCount = 0;
        for (int row = 0; row < observationDimension; ++row) {
            final double value = observationModel.getObservations().getParameterValue(row, timeIndex);
            if (!Double.isNaN(value)) {
                observedIndices[observedCount] = row;
                observedVector[observedCount] = value;
                ++observedCount;
            }
        }

        if (observedCount == 0) {
            return false;
        }

        copyObservedDesignRows();
        copyObservedNoiseBlock();

        final double noiseLogDet = invertObservedNoiseCovariance();
        buildObservationPrecisionContribution();
        buildObservationInformationContribution();
        potentialLogNormalizer = observationPotentialLogNormalizer(noiseLogDet);
        return true;
    }

    double applyTo(final CanonicalGaussianState predictedState,
                   final CanonicalGaussianState filteredState,
                   final CanonicalGaussianMessageOps.Workspace workspace) {
        addMatrices(predictedState.precision, observationPrecisionContribution, filteredState.precision);
        addVectors(predictedState.information, observationInformationContribution, filteredState.information);
        filteredState.logNormalizer = CanonicalGaussianMessageOps.normalizedLogNormalizer(filteredState, workspace);
        return filteredState.logNormalizer
                - predictedState.logNormalizer
                - potentialLogNormalizer;
    }

    private void copyObservedDesignRows() {
        for (int row = 0; row < observedCount; ++row) {
            final int sourceRow = observedIndices[row];
            System.arraycopy(fullDesignMatrix, sourceRow * stateDimension,
                    observedDesignMatrix, row * stateDimension,
                    stateDimension);
        }
    }

    private void copyObservedNoiseBlock() {
        for (int row = 0; row < observedCount; ++row) {
            final int sourceRow = observedIndices[row];
            final int targetRowOffset = row * observedCount;
            final int sourceRowOffset = sourceRow * observationDimension;
            for (int col = 0; col < observedCount; ++col) {
                observedNoiseCovariance[targetRowOffset + col] =
                        fullNoiseCovariance[sourceRowOffset + observedIndices[col]];
            }
        }
    }

    private double invertObservedNoiseCovariance() {
        if (!MatrixOps.tryCholesky(observedNoiseCovariance, flatCholeskyWorkspace, observedCount)) {
            throw new IllegalArgumentException("Observation noise submatrix is not positive definite");
        }
        return MatrixOps.invertFromCholesky(
                flatCholeskyWorkspace, flatLowerInverseWorkspace, observedNoisePrecision, observedCount);
    }

    private void buildObservationPrecisionContribution() {
        MatrixOps.matMul(observedNoisePrecision, observedDesignMatrix, noisePrecisionTimesDesign,
                observedCount, observedCount, stateDimension);
        MatrixOps.matMulTransposedLeft(observedDesignMatrix, noisePrecisionTimesDesign,
                observationPrecisionContribution, observedCount, stateDimension);
        MatrixOps.symmetrize(observationPrecisionContribution, stateDimension);
    }

    private void buildObservationInformationContribution() {
        MatrixOps.matVec(observedNoisePrecision, observedVector, noisePrecisionTimesObservation,
                observedCount, observedCount);
        for (int i = 0; i < stateDimension; ++i) {
            double sum = 0.0;
            for (int row = 0; row < observedCount; ++row) {
                sum += observedDesignMatrix[row * stateDimension + i] * noisePrecisionTimesObservation[row];
            }
            observationInformationContribution[i] = sum;
        }
    }

    private double observationPotentialLogNormalizer(final double logDetNoise) {
        final double quadratic = MatrixOps.quadraticForm(
                observedVector, observedNoisePrecision, observedCount, observationWorkspace);
        return 0.5 * (observedCount * MatrixOps.LOG_TWO_PI + logDetNoise + quadratic);
    }

    private static void addMatrices(final double[] left,
                                    final double[] right,
                                    final double[] out) {
        for (int i = 0; i < left.length; ++i) {
            out[i] = left[i] + right[i];
        }
    }

    private static void addVectors(final double[] left,
                                   final double[] right,
                                   final double[] out) {
        for (int i = 0; i < left.length; ++i) {
            out[i] = left[i] + right[i];
        }
    }
}
