package dr.inference.timeseries.engine.gaussian;

import dr.evomodel.treedatalikelihood.continuous.canonical.message.GaussianMatrixOps;

import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalBranchMessageContribution;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalBranchMessageContributionUtils;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalLocalTransitionAdjoints;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalTransitionAdjointUtils;

import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import dr.inference.timeseries.core.TimeGrid;
import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.continuous.ou.orthogonalblockdiagonal.OrthogonalBlockCanonicalParameterization;
import dr.evomodel.treedatalikelihood.continuous.canonical.math.GaussianFormConverter;
import dr.inference.timeseries.representation.GaussianTransitionRepresentation;

/**
 * Canonical-path analytical gradient of the log-likelihood with respect to the OU stationary mean.
 */
public final class CanonicalStationaryMeanGradientFormula implements CanonicalGradientFormula {

    private final Parameter stationaryMeanParameter;
    private final MatrixParameter initialCovarianceParameter;
    private final int stateDimension;
    private final OUProcessModel processModel;

    private final CanonicalBranchMessageContribution localContribution;
    private final CanonicalBranchMessageContributionUtils.Workspace contributionWorkspace;
    private final CanonicalLocalTransitionAdjoints localAdjoints;
    private final CanonicalTransitionAdjointUtils.Workspace transitionWorkspace;
    private final double[] currentMean;
    private final double[] smoothedInitialMean;
    private final double[] denseGradient;
    private final double[] stateDiff;
    private final double[][] initialCovInv;
    private final double[][] smoothedInitialCovariance;
    private final CanonicalOrthogonalBlockGradientCache orthogonalBlockGradientCache;

    public CanonicalStationaryMeanGradientFormula(final Parameter stationaryMeanParameter,
                                                  final MatrixParameter initialCovarianceParameter,
                                                  final int stateDimension) {
        this(null, stationaryMeanParameter, initialCovarianceParameter, stateDimension);
    }

    public CanonicalStationaryMeanGradientFormula(final OUProcessModel processModel,
                                                  final Parameter stationaryMeanParameter,
                                                  final MatrixParameter initialCovarianceParameter,
                                                  final int stateDimension) {
        this(processModel, stationaryMeanParameter, initialCovarianceParameter, stateDimension, null);
    }

    public CanonicalStationaryMeanGradientFormula(final OUProcessModel processModel,
                                                  final Parameter stationaryMeanParameter,
                                                  final MatrixParameter initialCovarianceParameter,
                                                  final int stateDimension,
                                                  final CanonicalOrthogonalBlockGradientCache orthogonalBlockGradientCache) {
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
        this.orthogonalBlockGradientCache = orthogonalBlockGradientCache;

        this.localContribution = new CanonicalBranchMessageContribution(stateDimension);
        this.contributionWorkspace = new CanonicalBranchMessageContributionUtils.Workspace(stateDimension);
        this.localAdjoints = new CanonicalLocalTransitionAdjoints(stateDimension);
        this.transitionWorkspace = new CanonicalTransitionAdjointUtils.Workspace(stateDimension);
        this.currentMean = new double[stateDimension];
        this.smoothedInitialMean = new double[stateDimension];
        this.denseGradient = new double[stateDimension];
        this.stateDiff = new double[stateDimension];
        this.initialCovInv = new double[stateDimension][stateDimension];
        this.smoothedInitialCovariance = new double[stateDimension][stateDimension];
    }

    @Override
    public boolean supportsParameter(final Parameter parameter) {
        return parameter == stationaryMeanParameter;
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
        if (stationaryMeanParameter.getDimension() == 1) {
            if (orthogonalBlockGradientCache != null
                    && orthogonalBlockGradientCache.supportsMeanParameter(parameter)
                    && branchGradientCache != null) {
                return orthogonalBlockGradientCache.getMeanGradient(
                        trajectory, branchGradientCache, timeGrid);
            }
            return computeScalarGradient(trajectory, branchGradientCache, timeGrid);
        }
        if (orthogonalBlockGradientCache != null
                && orthogonalBlockGradientCache.supportsMeanParameter(parameter)
                && branchGradientCache != null) {
            return orthogonalBlockGradientCache.getMeanGradient(
                    trajectory, branchGradientCache, timeGrid);
        }

        zero(denseGradient);
        fillCurrentMean(currentMean);

        final int timeCount = trajectory.timeCount;
        final OrthogonalBlockCanonicalParameterization orthogonalParameterization =
                getOrthogonalParameterizationIfAvailable();
        if (branchGradientCache != null) {
            branchGradientCache.ensure(trajectory);
        }
        for (int t = 0; t < timeCount - 1; ++t) {
            final CanonicalLocalTransitionAdjoints adjoints =
                    localAdjoints(t, trajectory, branchGradientCache);
            if (orthogonalParameterization != null) {
                orthogonalParameterization.accumulateMeanGradient(
                        timeGrid.getDelta(t, t + 1),
                        adjoints.dLogL_df,
                        denseGradient);
            } else {
                accumulateBranchMeanGradient(
                        branchGradientCache == null
                                ? transitionWorkspace.transitionMatrix
                                : branchGradientCache.getTransitionMatrix(t),
                        adjoints.dLogL_df,
                        denseGradient);
            }
        }

        GaussianFormConverter.fillMomentsFromState(
                trajectory.smoothedStates[0],
                smoothedInitialMean,
                smoothedInitialCovariance);
        GaussianMatrixOps.copyMatrix(initialCovarianceParameter.getParameterAsMatrix(), initialCovInv);
        final GaussianMatrixOps.CholeskyFactor initialChol =
                GaussianMatrixOps.cholesky(initialCovInv);
        GaussianMatrixOps.invertPositiveDefiniteFromCholesky(initialCovInv, initialChol);
        for (int i = 0; i < stateDimension; ++i) {
            stateDiff[i] = smoothedInitialMean[i] - currentMean[i];
        }
        final double[] initialGradient = localAdjoints.dLogL_df;
        GaussianMatrixOps.multiplyMatrixVector(initialCovInv, stateDiff, initialGradient);
        for (int i = 0; i < stateDimension; ++i) {
            denseGradient[i] += initialGradient[i];
        }

        return projectToParameterDimension(denseGradient);
    }

    private double[] computeScalarGradient(final CanonicalForwardTrajectory trajectory,
                                           final CanonicalBranchGradientCache branchGradientCache,
                                           final TimeGrid timeGrid) {
        double scalarGradient = 0.0;
        final double meanValue = stationaryMeanParameter.getParameterValue(0);
        final int timeCount = trajectory.timeCount;
        final OrthogonalBlockCanonicalParameterization orthogonalParameterization =
                getOrthogonalParameterizationIfAvailable();
        if (branchGradientCache != null) {
            branchGradientCache.ensure(trajectory);
        }

        for (int t = 0; t < timeCount - 1; ++t) {
            final CanonicalLocalTransitionAdjoints adjoints =
                    localAdjoints(t, trajectory, branchGradientCache);
            if (orthogonalParameterization != null) {
                scalarGradient += orthogonalParameterization.accumulateScalarMeanGradient(
                        timeGrid.getDelta(t, t + 1),
                        adjoints.dLogL_df);
            } else {
                scalarGradient += accumulateScalarBranchMeanGradient(
                        branchGradientCache == null
                                ? transitionWorkspace.transitionMatrix
                                : branchGradientCache.getTransitionMatrix(t),
                        adjoints.dLogL_df);
            }
        }

        GaussianFormConverter.fillMomentsFromState(
                trajectory.smoothedStates[0],
                smoothedInitialMean,
                smoothedInitialCovariance);
        GaussianMatrixOps.copyMatrix(initialCovarianceParameter.getParameterAsMatrix(), initialCovInv);
        final GaussianMatrixOps.CholeskyFactor initialChol =
                GaussianMatrixOps.cholesky(initialCovInv);
        GaussianMatrixOps.invertPositiveDefiniteFromCholesky(initialCovInv, initialChol);
        for (int i = 0; i < stateDimension; ++i) {
            stateDiff[i] = smoothedInitialMean[i] - meanValue;
        }
        final double[] initialGradient = localAdjoints.dLogL_df;
        GaussianMatrixOps.multiplyMatrixVector(initialCovInv, stateDiff, initialGradient);
        for (int i = 0; i < stateDimension; ++i) {
            scalarGradient += initialGradient[i];
        }

        return new double[]{scalarGradient};
    }

    private OrthogonalBlockCanonicalParameterization getOrthogonalParameterizationIfAvailable() {
        if (processModel == null) {
            return null;
        }
        if (processModel.getSelectionMatrixParameterization()
                instanceof OrthogonalBlockCanonicalParameterization) {
            return (OrthogonalBlockCanonicalParameterization)
                    processModel.getSelectionMatrixParameterization();
        }
        return null;
    }

    @Override
    public void makeDirty() {
        if (orthogonalBlockGradientCache != null) {
            orthogonalBlockGradientCache.makeDirty();
        }
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

    private static void accumulateBranchMeanGradient(final double[] transitionMatrix,
                                                     final double[] dLogL_df,
                                                     final double[] accumulator) {
        final int d = dLogL_df.length;
        for (int j = 0; j < d; ++j) {
            double sum = dLogL_df[j];
            for (int i = 0; i < d; ++i) {
                sum -= transitionMatrix[i * d + j] * dLogL_df[i];
            }
            accumulator[j] += sum;
        }
    }

    private static double accumulateScalarBranchMeanGradient(final double[] transitionMatrix,
                                                             final double[] dLogL_df) {
        final int d = dLogL_df.length;
        double accumulator = 0.0;
        for (int j = 0; j < d; ++j) {
            double sum = dLogL_df[j];
            for (int i = 0; i < d; ++i) {
                sum -= transitionMatrix[i * d + j] * dLogL_df[i];
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

    private CanonicalLocalTransitionAdjoints localAdjoints(final int branchIndex,
                                                           final CanonicalForwardTrajectory trajectory,
                                                           final CanonicalBranchGradientCache branchGradientCache) {
        if (branchGradientCache != null) {
            return branchGradientCache.getAdjoints(branchIndex);
        }
        CanonicalBranchMessageContributionUtils.fillFromPairState(
                trajectory.branchPairStates[branchIndex],
                contributionWorkspace,
                localContribution);
        CanonicalTransitionAdjointUtils.fillFromCanonicalTransition(
                trajectory.transitions[branchIndex],
                localContribution,
                transitionWorkspace,
                localAdjoints);
        return localAdjoints;
    }
}
