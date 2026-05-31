package dr.inference.timeseries.engine.kalman.formula;

import dr.evomodel.treedatalikelihood.continuous.canonical.message.GaussianMatrixOps;

import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalLocalTransitionAdjoints;

import dr.inference.model.Parameter;
import dr.inference.timeseries.core.TimeGrid;
import dr.evomodel.continuous.ou.DiffusionMatrixParameterization;
import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.continuous.ou.blockdiagonal.BlockDiagonalNativeCanonicalParameterization;
import dr.inference.timeseries.engine.kalman.CanonicalBranchGradientCache;
import dr.inference.timeseries.engine.kalman.CanonicalForwardTrajectory;
import dr.inference.timeseries.representation.GaussianTransitionRepresentation;

/**
 * Diffusion gradient formula driven by branch-local canonical posterior states.
 */
public final class CanonicalDiffusionMatrixGradientFormula implements CanonicalGradientFormula {

    private final DiffusionMatrixParameterization diffusionParameterization;
    private final int stateDimension;
    private final OUProcessModel processModel;

    private final CanonicalBranchAdjointProvider branchAdjoints;
    private final double[][] covarianceAdjointScratch;
    private final CanonicalBlockDiagonalGradientCache blockDiagonalGradientCache;

    public CanonicalDiffusionMatrixGradientFormula(final DiffusionMatrixParameterization diffusionParameterization,
                                                   final int stateDimension) {
        this(null, diffusionParameterization, stateDimension);
    }

    public CanonicalDiffusionMatrixGradientFormula(final OUProcessModel processModel,
                                                   final DiffusionMatrixParameterization diffusionParameterization,
                                                   final int stateDimension) {
        this(processModel, diffusionParameterization, stateDimension, null);
    }

    public CanonicalDiffusionMatrixGradientFormula(final OUProcessModel processModel,
                                                   final DiffusionMatrixParameterization diffusionParameterization,
                                                   final int stateDimension,
                                                   final CanonicalBlockDiagonalGradientCache blockDiagonalGradientCache) {
        if (diffusionParameterization == null) {
            throw new IllegalArgumentException("diffusionParameterization must not be null");
        }
        if (stateDimension < 1) {
            throw new IllegalArgumentException("stateDimension must be at least 1");
        }
        this.processModel = processModel;
        this.diffusionParameterization = diffusionParameterization;
        this.stateDimension = stateDimension;
        this.blockDiagonalGradientCache = blockDiagonalGradientCache;

        this.branchAdjoints = new CanonicalBranchAdjointProvider(stateDimension);
        this.covarianceAdjointScratch = new double[stateDimension][stateDimension];
    }

    @Override
    public boolean supportsParameter(final Parameter parameter) {
        return diffusionParameterization.supportsParameter(parameter);
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
        final int d = stateDimension;
        final int T = trajectory.timeCount;
        if (blockDiagonalGradientCache != null
                && blockDiagonalGradientCache.supportsDiffusionParameter(parameter)
                && branchGradientCache != null) {
            return blockDiagonalGradientCache.getDiffusionGradient(
                    parameter, trajectory, branchGradientCache, repr, timeGrid);
        }
        final double[] gradientAccumulator = new double[d * d];
        final BlockDiagonalNativeCanonicalParameterization blockParameterization =
                BlockDiagonalFormulaSupport.nativeParameterization(processModel);
        branchAdjoints.ensure(trajectory, branchGradientCache);

        for (int t = 0; t < T - 1; ++t) {
            final CanonicalLocalTransitionAdjoints adjoints =
                    branchAdjoints.localAdjoints(t, trajectory, branchGradientCache);
            if (blockParameterization != null) {
                blockParameterization.accumulateDiffusionGradientFlat(
                        processModel.getDiffusionMatrix(),
                        timeGrid.getDelta(t, t + 1),
                        adjoints.dLogL_dOmega,
                        false,
                        gradientAccumulator);
            } else {
                GaussianMatrixOps.copyFlatToMatrix(adjoints.dLogL_dOmega, covarianceAdjointScratch, d);
                repr.accumulateDiffusionGradient(
                        t, t + 1, timeGrid, covarianceAdjointScratch, gradientAccumulator);
            }
        }

        return diffusionParameterization.pullBackGradient(parameter, gradientAccumulator);
    }

    @Override
    public void makeDirty() {
        if (blockDiagonalGradientCache != null) {
            blockDiagonalGradientCache.makeDirty();
        }
    }

}
