package dr.inference.timeseries.engine.gaussian;

import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import dr.inference.timeseries.core.TimeGrid;
import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.continuous.ou.OrthogonalBlockDiagonalSelectionMatrixParameterization;
import dr.inference.timeseries.representation.GaussianTransitionRepresentation;

/**
 * Analytical gradient of the log-likelihood with respect to the OU stationary/optimal mean.
 *
 * <p>This parameter affects:
 * <ul>
 *   <li>the branch offset {@code f_t = (I - F_t) mu}</li>
 *   <li>the initial latent mean {@code x_0 ~ N(mu, P0)}</li>
 * </ul>
 * so the total gradient is the sum of branch contributions and the initial-state term.
 */
public final class StationaryMeanGradientFormula implements GradientFormula {

    private final Parameter stationaryMeanParameter;
    private final MatrixParameter initialCovarianceParameter;
    private final int stateDimension;
    private final OUProcessModel processModel;

    private final double[] meanResidual;
    private final double[][] stepCovInv;
    private final double[] dLogL_df;
    private final double[] currentMean;
    private final double[] stateDiff;
    private final double[] denseGradient;

    public StationaryMeanGradientFormula(final Parameter stationaryMeanParameter,
                                         final MatrixParameter initialCovarianceParameter,
                                         final int stateDimension) {
        this(null, stationaryMeanParameter, initialCovarianceParameter, stateDimension);
    }

    public StationaryMeanGradientFormula(final OUProcessModel processModel,
                                         final Parameter stationaryMeanParameter,
                                         final MatrixParameter initialCovarianceParameter,
                                         final int stateDimension) {
        if (stationaryMeanParameter == null) {
            throw new IllegalArgumentException("stationaryMeanParameter must not be null");
        }
        if (initialCovarianceParameter == null) {
            throw new IllegalArgumentException("initialCovarianceParameter must not be null");
        }
        if (stateDimension < 1) {
            throw new IllegalArgumentException("stateDimension must be at least 1");
        }
        this.processModel = processModel;
        this.stationaryMeanParameter = stationaryMeanParameter;
        this.initialCovarianceParameter = initialCovarianceParameter;
        this.stateDimension = stateDimension;

        this.meanResidual = new double[stateDimension];
        this.stepCovInv = new double[stateDimension][stateDimension];
        this.dLogL_df = new double[stateDimension];
        this.currentMean = new double[stateDimension];
        this.stateDiff = new double[stateDimension];
        this.denseGradient = new double[stateDimension];
    }

    @Override
    public boolean supportsParameter(final Parameter parameter) {
        return parameter == stationaryMeanParameter;
    }

    @Override
    public double[] computeGradient(final Parameter parameter,
                                    final BranchSmootherStats[] smootherStats,
                                    final ForwardTrajectory trajectory,
                                    final GaussianTransitionRepresentation repr,
                                    final TimeGrid timeGrid) {
        if (stationaryMeanParameter.getDimension() == 1) {
            return computeScalarGradient(smootherStats, trajectory, timeGrid);
        }

        zero(denseGradient);
        fillCurrentMean(currentMean);

        final int timeCount = trajectory.timeCount;
        final OrthogonalBlockDiagonalSelectionMatrixParameterization orthogonalParameterization =
                getOrthogonalParameterizationIfAvailable();
        for (int t = 0; t < timeCount - 1; ++t) {
            final BranchSmootherStats curr = smootherStats[t];
            final BranchSmootherStats next = smootherStats[t + 1];
            final double[][] transitionMatrix = trajectory.transitionMatrices[t];
            final double[] transitionOffset = trajectory.transitionOffsets[t];

            KalmanLikelihoodEngine.multiplyMatrixVector(transitionMatrix, curr.smoothedMean, meanResidual);
            for (int i = 0; i < stateDimension; ++i) {
                meanResidual[i] = next.smoothedMean[i] - meanResidual[i] - transitionOffset[i];
            }

            KalmanLikelihoodEngine.copyMatrix(trajectory.stepCovariances[t], stepCovInv);
            final KalmanLikelihoodEngine.CholeskyFactor stepChol =
                    KalmanLikelihoodEngine.cholesky(stepCovInv);
            KalmanLikelihoodEngine.invertPositiveDefiniteFromCholesky(stepCovInv, stepChol);
            KalmanLikelihoodEngine.multiplyMatrixVector(stepCovInv, meanResidual, dLogL_df);

            if (orthogonalParameterization != null) {
                orthogonalParameterization.accumulateMeanGradient(
                        timeGrid.getDelta(t, t + 1),
                        dLogL_df,
                        denseGradient);
            } else {
                accumulateBranchMeanGradient(transitionMatrix, dLogL_df, denseGradient);
            }
        }

        KalmanLikelihoodEngine.copyMatrix(initialCovarianceParameter.getParameterAsMatrix(), stepCovInv);
        final KalmanLikelihoodEngine.CholeskyFactor initialChol =
                KalmanLikelihoodEngine.cholesky(stepCovInv);
        KalmanLikelihoodEngine.invertPositiveDefiniteFromCholesky(stepCovInv, initialChol);
        for (int i = 0; i < stateDimension; ++i) {
            stateDiff[i] = smootherStats[0].smoothedMean[i] - currentMean[i];
        }
        KalmanLikelihoodEngine.multiplyMatrixVector(stepCovInv, stateDiff, dLogL_df);
        for (int i = 0; i < stateDimension; ++i) {
            denseGradient[i] += dLogL_df[i];
        }

        return projectToParameterDimension(denseGradient);
    }

    private double[] computeScalarGradient(final BranchSmootherStats[] smootherStats,
                                           final ForwardTrajectory trajectory,
                                           final TimeGrid timeGrid) {
        double scalarGradient = 0.0;
        final double meanValue = stationaryMeanParameter.getParameterValue(0);
        final int timeCount = trajectory.timeCount;
        final OrthogonalBlockDiagonalSelectionMatrixParameterization orthogonalParameterization =
                getOrthogonalParameterizationIfAvailable();

        for (int t = 0; t < timeCount - 1; ++t) {
            final BranchSmootherStats curr = smootherStats[t];
            final BranchSmootherStats next = smootherStats[t + 1];
            final double[][] transitionMatrix = trajectory.transitionMatrices[t];
            final double[] transitionOffset = trajectory.transitionOffsets[t];

            KalmanLikelihoodEngine.multiplyMatrixVector(transitionMatrix, curr.smoothedMean, meanResidual);
            for (int i = 0; i < stateDimension; ++i) {
                meanResidual[i] = next.smoothedMean[i] - meanResidual[i] - transitionOffset[i];
            }

            KalmanLikelihoodEngine.copyMatrix(trajectory.stepCovariances[t], stepCovInv);
            final KalmanLikelihoodEngine.CholeskyFactor stepChol =
                    KalmanLikelihoodEngine.cholesky(stepCovInv);
            KalmanLikelihoodEngine.invertPositiveDefiniteFromCholesky(stepCovInv, stepChol);
            KalmanLikelihoodEngine.multiplyMatrixVector(stepCovInv, meanResidual, dLogL_df);

            if (orthogonalParameterization != null) {
                scalarGradient += orthogonalParameterization.accumulateScalarMeanGradient(
                        timeGrid.getDelta(t, t + 1),
                        dLogL_df);
            } else {
                scalarGradient += accumulateScalarBranchMeanGradient(transitionMatrix, dLogL_df);
            }
        }

        KalmanLikelihoodEngine.copyMatrix(initialCovarianceParameter.getParameterAsMatrix(), stepCovInv);
        final KalmanLikelihoodEngine.CholeskyFactor initialChol =
                KalmanLikelihoodEngine.cholesky(stepCovInv);
        KalmanLikelihoodEngine.invertPositiveDefiniteFromCholesky(stepCovInv, initialChol);
        for (int i = 0; i < stateDimension; ++i) {
            stateDiff[i] = smootherStats[0].smoothedMean[i] - meanValue;
        }
        KalmanLikelihoodEngine.multiplyMatrixVector(stepCovInv, stateDiff, dLogL_df);
        for (int i = 0; i < stateDimension; ++i) {
            scalarGradient += dLogL_df[i];
        }

        return new double[]{scalarGradient};
    }

    private OrthogonalBlockDiagonalSelectionMatrixParameterization getOrthogonalParameterizationIfAvailable() {
        if (processModel == null) {
            return null;
        }
        if (processModel.getSelectionMatrixParameterization()
                instanceof OrthogonalBlockDiagonalSelectionMatrixParameterization) {
            return (OrthogonalBlockDiagonalSelectionMatrixParameterization)
                    processModel.getSelectionMatrixParameterization();
        }
        return null;
    }

    private void fillCurrentMean(final double[] out) {
        if (stationaryMeanParameter.getDimension() == 1) {
            final double value = stationaryMeanParameter.getParameterValue(0);
            for (int i = 0; i < stateDimension; ++i) {
                out[i] = value;
            }
        } else if (stationaryMeanParameter.getDimension() == stateDimension) {
            for (int i = 0; i < stateDimension; ++i) {
                out[i] = stationaryMeanParameter.getParameterValue(i);
            }
        } else {
            throw new IllegalStateException("stationaryMean dimension must be 1 or stateDimension");
        }
    }

    private double[] projectToParameterDimension(final double[] denseSource) {
        if (stationaryMeanParameter.getDimension() == stateDimension) {
            return denseSource.clone();
        }
        double sum = 0.0;
        for (int i = 0; i < stateDimension; ++i) {
            sum += denseSource[i];
        }
        return new double[]{sum};
    }

    private static void accumulateBranchMeanGradient(final double[][] transitionMatrix,
                                                     final double[] dLogL_df,
                                                     final double[] accumulator) {
        final int d = dLogL_df.length;
        for (int j = 0; j < d; ++j) {
            double sum = dLogL_df[j];
            for (int i = 0; i < d; ++i) {
                sum -= transitionMatrix[i][j] * dLogL_df[i];
            }
            accumulator[j] += sum;
        }
    }

    private static double accumulateScalarBranchMeanGradient(final double[][] transitionMatrix,
                                                             final double[] dLogL_df) {
        final int d = dLogL_df.length;
        double accumulator = 0.0;
        for (int j = 0; j < d; ++j) {
            double sum = dLogL_df[j];
            for (int i = 0; i < d; ++i) {
                sum -= transitionMatrix[i][j] * dLogL_df[i];
            }
            accumulator += sum;
        }
        return accumulator;
    }

    private static void zero(final double[] vector) {
        for (int i = 0; i < vector.length; ++i) {
            vector[i] = 0.0;
        }
    }
}
