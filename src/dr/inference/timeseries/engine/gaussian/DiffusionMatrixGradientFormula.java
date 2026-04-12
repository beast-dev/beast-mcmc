package dr.inference.timeseries.engine.gaussian;

import dr.inference.model.Parameter;
import dr.inference.timeseries.core.TimeGrid;
import dr.inference.timeseries.gaussian.DiffusionMatrixParameterization;
import dr.inference.timeseries.representation.GaussianTransitionRepresentation;

/**
 * Analytical gradient of the Kalman log-likelihood with respect to the diffusion
 * matrix Q (or native parameters that determine Q).
 */
public final class DiffusionMatrixGradientFormula implements GradientFormula {

    private final DiffusionMatrixParameterization diffusionParameterization;
    private final int stateDimension;

    private final double[] meanResidual;
    private final double[][] crossCov;
    private final double[][] stepCovInv;
    private final double[][] tempDxD;
    private final double[][] tempDxD2;
    private final double[][] residualSecondMoment;
    private final double[][] dLogL_dV;

    public DiffusionMatrixGradientFormula(final DiffusionMatrixParameterization diffusionParameterization,
                                          final int stateDimension) {
        if (diffusionParameterization == null) {
            throw new IllegalArgumentException("diffusionParameterization must not be null");
        }
        if (stateDimension < 1) {
            throw new IllegalArgumentException("stateDimension must be at least 1");
        }
        this.diffusionParameterization = diffusionParameterization;
        this.stateDimension = stateDimension;

        meanResidual = new double[stateDimension];
        crossCov = new double[stateDimension][stateDimension];
        stepCovInv = new double[stateDimension][stateDimension];
        tempDxD = new double[stateDimension][stateDimension];
        tempDxD2 = new double[stateDimension][stateDimension];
        residualSecondMoment = new double[stateDimension][stateDimension];
        dLogL_dV = new double[stateDimension][stateDimension];
    }

    @Override
    public boolean supportsParameter(final Parameter parameter) {
        return diffusionParameterization.supportsParameter(parameter);
    }

    @Override
    public double[] computeGradient(final Parameter parameter,
                                    final BranchSmootherStats[] smootherStats,
                                    final ForwardTrajectory trajectory,
                                    final GaussianTransitionRepresentation repr,
                                    final TimeGrid timeGrid) {
        final int T = trajectory.timeCount;
        final int d = stateDimension;
        final double[] gradientAccumulator = new double[d * d];

        for (int t = 0; t < T - 1; ++t) {
            final BranchSmootherStats curr = smootherStats[t];
            final BranchSmootherStats next = smootherStats[t + 1];
            final double[][] F_t = trajectory.transitionMatrices[t];
            final double[] f_t = trajectory.transitionOffsets[t];

            KalmanLikelihoodEngine.multiplyMatrixMatrixTransposedRight(
                    next.smoothedCovariance, curr.smootherGain, crossCov);

            KalmanLikelihoodEngine.multiplyMatrixVector(F_t, curr.smoothedMean, meanResidual);
            for (int i = 0; i < d; ++i) {
                meanResidual[i] = next.smoothedMean[i] - meanResidual[i] - f_t[i];
            }

            KalmanLikelihoodEngine.copyMatrix(trajectory.stepCovariances[t], stepCovInv);
            final KalmanLikelihoodEngine.CholeskyFactor stepChol =
                    KalmanLikelihoodEngine.cholesky(stepCovInv);
            KalmanLikelihoodEngine.invertPositiveDefiniteFromCholesky(stepCovInv, stepChol);

            KalmanLikelihoodEngine.copyMatrix(next.smoothedCovariance, residualSecondMoment);
            KalmanLikelihoodEngine.multiplyMatrixMatrixTransposedRight(F_t, crossCov, tempDxD);
            for (int i = 0; i < d; ++i) {
                for (int j = 0; j < d; ++j) {
                    residualSecondMoment[i][j] -= tempDxD[i][j];
                }
            }
            KalmanLikelihoodEngine.multiplyMatrixMatrixTransposedRight(crossCov, F_t, tempDxD);
            for (int i = 0; i < d; ++i) {
                for (int j = 0; j < d; ++j) {
                    residualSecondMoment[i][j] -= tempDxD[i][j];
                }
            }
            KalmanLikelihoodEngine.multiplyMatrixMatrix(F_t, curr.smoothedCovariance, tempDxD);
            KalmanLikelihoodEngine.multiplyMatrixMatrixTransposedRight(tempDxD, F_t, tempDxD2);
            for (int i = 0; i < d; ++i) {
                for (int j = 0; j < d; ++j) {
                    residualSecondMoment[i][j] += tempDxD2[i][j];
                    residualSecondMoment[i][j] += meanResidual[i] * meanResidual[j];
                }
            }

            KalmanLikelihoodEngine.multiplyMatrixMatrix(stepCovInv, residualSecondMoment, tempDxD);
            KalmanLikelihoodEngine.multiplyMatrixMatrix(tempDxD, stepCovInv, dLogL_dV);
            for (int i = 0; i < d; ++i) {
                for (int j = 0; j < d; ++j) {
                    dLogL_dV[i][j] = 0.5 * (dLogL_dV[i][j] - stepCovInv[i][j]);
                }
            }

            repr.accumulateDiffusionGradient(t, t + 1, timeGrid, dLogL_dV, gradientAccumulator);
        }

        if (parameter == diffusionParameterization.getMatrixParameter()) {
            return gradientAccumulator;
        }
        return diffusionParameterization.pullBackGradient(parameter, gradientAccumulator);
    }
}
