package dr.inference.timeseries.engine.kalman.formula;

import dr.evomodel.treedatalikelihood.continuous.canonical.math.MatrixOps;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianState;
import dr.inference.model.DiagonalMatrix;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import dr.inference.model.ScaledIdentityMatrixParameter;
import dr.inference.timeseries.core.TimeGrid;
import dr.inference.timeseries.engine.kalman.CanonicalBranchGradientCache;
import dr.inference.timeseries.engine.kalman.CanonicalForwardTrajectory;
import dr.inference.timeseries.engine.kalman.MomentBranchSmootherStats;
import dr.inference.timeseries.engine.kalman.MomentForwardTrajectory;
import dr.inference.timeseries.model.gaussian.LinearGaussianObservationModel;
import dr.inference.timeseries.representation.GaussianTransitionRepresentation;

/**
 * Analytical gradient of the Kalman log likelihood with respect to Gaussian
 * observation-noise variances for diagonal or scalar-identity covariance.
 */
public final class DiagonalObservationNoiseGradientFormula
        implements MomentGradientFormula, CanonicalGradientFormula {

    private enum NoiseMode {
        DIAGONAL,
        SCALAR_IDENTITY,
        UNSUPPORTED
    }

    private final LinearGaussianObservationModel observationModel;
    private final MatrixParameter noiseCovariance;
    private final Parameter gradientParameter;
    private final NoiseMode noiseMode;
    private final int stateDimension;
    private final int observationDimension;

    private final double[] designMatrix;
    private final double[] canonicalCovariance;
    private final double[] canonicalCholesky;
    private final double[] canonicalLowerInverse;
    private final double[] canonicalMean;

    public DiagonalObservationNoiseGradientFormula(final LinearGaussianObservationModel observationModel,
                                                   final int stateDimension) {
        if (observationModel == null) {
            throw new IllegalArgumentException("observationModel must not be null");
        }
        if (stateDimension < 1) {
            throw new IllegalArgumentException("stateDimension must be at least 1");
        }
        this.observationModel = observationModel;
        this.stateDimension = stateDimension;
        this.observationDimension = observationModel.getObservationDimension();

        final MatrixParameter noise = observationModel.getNoiseCovariance();
        if (noise instanceof ScaledIdentityMatrixParameter
                && ((ScaledIdentityMatrixParameter) noise).getRowDimension() == observationDimension) {
            this.noiseCovariance = noise;
            this.gradientParameter = ((ScaledIdentityMatrixParameter) noise).getScaleParameter();
            this.noiseMode = NoiseMode.SCALAR_IDENTITY;
        } else if (noise instanceof DiagonalMatrix
                && ((DiagonalMatrix) noise).getDiagonalParameter().getDimension() == observationDimension) {
            this.noiseCovariance = (DiagonalMatrix) noise;
            this.gradientParameter = ((DiagonalMatrix) noise).getDiagonalParameter();
            this.noiseMode = NoiseMode.DIAGONAL;
        } else {
            this.noiseCovariance = null;
            this.gradientParameter = null;
            this.noiseMode = NoiseMode.UNSUPPORTED;
        }

        this.designMatrix = new double[observationDimension * stateDimension];
        this.canonicalCovariance = new double[stateDimension * stateDimension];
        this.canonicalCholesky = new double[stateDimension * stateDimension];
        this.canonicalLowerInverse = new double[stateDimension * stateDimension];
        this.canonicalMean = new double[stateDimension];
    }

    @Override
    public boolean supportsParameter(final Parameter parameter) {
        return gradientParameter != null && parameter == gradientParameter;
    }

    @Override
    public double[] computeGradient(final Parameter parameter,
                                    final MomentBranchSmootherStats[] smootherStats,
                                    final MomentForwardTrajectory trajectory,
                                    final GaussianTransitionRepresentation repr,
                                    final TimeGrid timeGrid) {
        checkSupported(parameter);
        observationModel.fillDesignMatrixFlat(designMatrix, stateDimension);

        final double[] gradient = new double[gradientParameter.getDimension()];
        for (int timeIndex = 0; timeIndex < trajectory.timeCount; ++timeIndex) {
            accumulateTimeGradient(
                    timeIndex,
                    smootherStats[timeIndex].smoothedMean,
                    smootherStats[timeIndex].smoothedCovariance,
                    gradient);
        }
        return gradient;
    }

    @Override
    public double[] computeGradient(final Parameter parameter,
                                    final CanonicalForwardTrajectory trajectory,
                                    final GaussianTransitionRepresentation repr,
                                    final TimeGrid timeGrid) {
        return computeGradient(parameter, trajectory, null, repr, timeGrid);
    }

    @Override
    public double[] computeGradient(final Parameter parameter,
                                    final CanonicalForwardTrajectory trajectory,
                                    final CanonicalBranchGradientCache branchGradientCache,
                                    final GaussianTransitionRepresentation repr,
                                    final TimeGrid timeGrid) {
        checkSupported(parameter);
        observationModel.fillDesignMatrixFlat(designMatrix, stateDimension);

        final double[] gradient = new double[gradientParameter.getDimension()];
        for (int timeIndex = 0; timeIndex < trajectory.timeCount; ++timeIndex) {
            fillCanonicalMoments(trajectory.smoothedStates[timeIndex]);
            accumulateTimeGradient(timeIndex, canonicalMean, canonicalCovariance, gradient);
        }
        return gradient;
    }

    private void accumulateTimeGradient(final int timeIndex,
                                        final double[] smoothedMean,
                                        final double[] smoothedCovariance,
                                        final double[] gradient) {
        for (int row = 0; row < observationDimension; ++row) {
            final double observation = observationModel.getObservations().getParameterValue(row, timeIndex);
            if (Double.isNaN(observation)) {
                continue;
            }

            final double fitted = projectedMean(row, smoothedMean);
            final double residual = observation - fitted;
            final double expectedSquaredResidual =
                    residual * residual + projectedVariance(row, smoothedCovariance);
            final double variance = noiseCovariance.getParameterValue(row, row);
            if (!(variance > 0.0) || Double.isNaN(variance)) {
                throw new IllegalArgumentException("Observation noise variance must be positive");
            }
            final double contribution = -0.5 * (1.0 / variance
                    - expectedSquaredResidual / (variance * variance));
            if (noiseMode == NoiseMode.SCALAR_IDENTITY) {
                gradient[0] += contribution;
            } else {
                gradient[row] += contribution;
            }
        }
    }

    private double projectedMean(final int observationRow,
                                 final double[] smoothedMean) {
        double value = 0.0;
        final int designOffset = observationRow * stateDimension;
        for (int state = 0; state < stateDimension; ++state) {
            value += designMatrix[designOffset + state] * smoothedMean[state];
        }
        return value;
    }

    private double projectedVariance(final int observationRow,
                                     final double[] smoothedCovariance) {
        double value = 0.0;
        final int designOffset = observationRow * stateDimension;
        for (int left = 0; left < stateDimension; ++left) {
            final double hLeft = designMatrix[designOffset + left];
            if (hLeft == 0.0) {
                continue;
            }
            final int covarianceOffset = left * stateDimension;
            for (int right = 0; right < stateDimension; ++right) {
                value += hLeft
                        * smoothedCovariance[covarianceOffset + right]
                        * designMatrix[designOffset + right];
            }
        }
        return value;
    }

    private void fillCanonicalMoments(final CanonicalGaussianState smoothedState) {
        if (!MatrixOps.tryCholesky(smoothedState.precision, canonicalCholesky, stateDimension)) {
            throw new IllegalArgumentException("Smoothed canonical precision is not positive definite");
        }
        MatrixOps.invertFromCholesky(
                canonicalCholesky,
                canonicalLowerInverse,
                canonicalCovariance,
                stateDimension);
        MatrixOps.matVec(canonicalCovariance, smoothedState.information, canonicalMean, stateDimension);
    }

    private void checkSupported(final Parameter parameter) {
        if (!supportsParameter(parameter)) {
            throw new IllegalArgumentException(
                    "Observation-noise gradient supports only diagonal observation covariance "
                            + "parameters or scalar scaled-identity observation covariance parameters");
        }
    }
}
