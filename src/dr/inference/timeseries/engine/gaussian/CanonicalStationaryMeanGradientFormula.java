package dr.inference.timeseries.engine.gaussian;

import dr.evomodel.treedatalikelihood.continuous.canonical.message.GaussianMatrixOps;

import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalLocalTransitionAdjoints;

import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import dr.inference.timeseries.core.TimeGrid;
import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.continuous.ou.blockdiagonal.BlockDiagonalCanonicalParameterization;
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

    private final CanonicalBranchAdjointProvider branchAdjoints;
    private final double[] currentMean;
    private final double[] smoothedInitialMean;
    private final double[] denseGradient;
    private final double[] stateDiff;
    private final double[] initialMeanAdjoint;
    private final double[][] initialCovInv;
    private final double[][] smoothedInitialCovariance;
    private final CanonicalBlockDiagonalGradientCache blockDiagonalGradientCache;

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
                                                  final CanonicalBlockDiagonalGradientCache blockDiagonalGradientCache) {
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
        this.blockDiagonalGradientCache = blockDiagonalGradientCache;

        this.branchAdjoints = new CanonicalBranchAdjointProvider(stateDimension);
        this.currentMean = new double[stateDimension];
        this.smoothedInitialMean = new double[stateDimension];
        this.denseGradient = new double[stateDimension];
        this.stateDiff = new double[stateDimension];
        this.initialMeanAdjoint = new double[stateDimension];
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
            if (blockDiagonalGradientCache != null
                    && blockDiagonalGradientCache.supportsMeanParameter(parameter)
                    && branchGradientCache != null) {
                return blockDiagonalGradientCache.getMeanGradient(
                        trajectory, branchGradientCache, repr, timeGrid);
            }
            return computeScalarGradient(trajectory, branchGradientCache, timeGrid);
        }
        if (blockDiagonalGradientCache != null
                && blockDiagonalGradientCache.supportsMeanParameter(parameter)
                && branchGradientCache != null) {
            return blockDiagonalGradientCache.getMeanGradient(
                    trajectory, branchGradientCache, repr, timeGrid);
        }

        BlockDiagonalFormulaSupport.zero(denseGradient);
        BlockDiagonalFormulaSupport.fillCurrentMean(
                stationaryMeanParameter,
                stateDimension,
                currentMean);

        final int timeCount = trajectory.timeCount;
        final BlockDiagonalCanonicalParameterization blockParameterization =
                BlockDiagonalFormulaSupport.canonicalParameterization(processModel);
        branchAdjoints.ensure(trajectory, branchGradientCache);
        for (int t = 0; t < timeCount - 1; ++t) {
            final CanonicalLocalTransitionAdjoints adjoints =
                    branchAdjoints.localAdjoints(t, trajectory, branchGradientCache);
            if (blockParameterization != null) {
                blockParameterization.accumulateMeanGradient(
                        timeGrid.getDelta(t, t + 1),
                        adjoints.dLogL_df,
                        denseGradient);
            } else {
                BlockDiagonalFormulaSupport.accumulateBranchMeanGradientFlat(
                        branchAdjoints.transitionMatrix(t, branchGradientCache),
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
        GaussianMatrixOps.multiplyMatrixVector(initialCovInv, stateDiff, initialMeanAdjoint);
        for (int i = 0; i < stateDimension; ++i) {
            denseGradient[i] += initialMeanAdjoint[i];
        }

        return BlockDiagonalFormulaSupport.projectMeanGradient(
                stationaryMeanParameter,
                stateDimension,
                denseGradient);
    }

    private double[] computeScalarGradient(final CanonicalForwardTrajectory trajectory,
                                           final CanonicalBranchGradientCache branchGradientCache,
                                           final TimeGrid timeGrid) {
        double scalarGradient = 0.0;
        final double meanValue = stationaryMeanParameter.getParameterValue(0);
        final int timeCount = trajectory.timeCount;
        final BlockDiagonalCanonicalParameterization blockParameterization =
                BlockDiagonalFormulaSupport.canonicalParameterization(processModel);
        branchAdjoints.ensure(trajectory, branchGradientCache);

        for (int t = 0; t < timeCount - 1; ++t) {
            final CanonicalLocalTransitionAdjoints adjoints =
                    branchAdjoints.localAdjoints(t, trajectory, branchGradientCache);
            if (blockParameterization != null) {
                scalarGradient += blockParameterization.accumulateScalarMeanGradient(
                        timeGrid.getDelta(t, t + 1),
                        adjoints.dLogL_df);
            } else {
                scalarGradient += BlockDiagonalFormulaSupport.accumulateScalarBranchMeanGradientFlat(
                        branchAdjoints.transitionMatrix(t, branchGradientCache),
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
        GaussianMatrixOps.multiplyMatrixVector(initialCovInv, stateDiff, initialMeanAdjoint);
        for (int i = 0; i < stateDimension; ++i) {
            scalarGradient += initialMeanAdjoint[i];
        }

        return new double[]{scalarGradient};
    }
    @Override
    public void makeDirty() {
        if (blockDiagonalGradientCache != null) {
            blockDiagonalGradientCache.makeDirty();
        }
    }
}
