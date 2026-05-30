package dr.inference.timeseries.engine.gaussian;

import dr.evomodel.treedatalikelihood.continuous.canonical.message.GaussianMatrixOps;

import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalBranchMessageContribution;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalBranchMessageContributionUtils;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalLocalTransitionAdjoints;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalTransitionAdjointUtils;
import dr.evomodel.treedatalikelihood.continuous.canonical.gradient.CanonicalSelectionGradientProjector;

import dr.inference.model.AbstractBlockDiagonalTwoByTwoMatrixParameter;
import dr.inference.model.OrthogonalMatrixProvider;
import dr.inference.model.Parameter;
import dr.inference.timeseries.core.TimeGrid;
import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.continuous.ou.orthogonalblockdiagonal.BlockDiagonalNativeCanonicalParameterization;
import dr.inference.timeseries.representation.GaussianTransitionRepresentation;

import java.util.Arrays;

/**
 * Selection-matrix gradient formula driven by branch-local canonical posterior states.
 *
 * <p>This implementation is the first step toward native canonical branch adjoints:
 * for each branch it constructs the joint posterior over {@code (x_t, x_{t+1})} in
 * canonical form, then derives the branch moments needed to build the usual
 * transition-space sensitivities {@code dL/dF}, {@code dL/df}, and {@code dL/dV}.
 */
public final class CanonicalSelectionMatrixGradientFormula implements CanonicalGradientFormula {

    private static final String DISABLE_ORTHOGONAL_NATIVE_SELECTION_PATH_PROPERTY =
            "beast.experimental.disableOrthogonalNativeSelectionPath";
    private static final String DISABLE_BLOCK_DIAGONAL_NATIVE_SELECTION_PATH_PROPERTY =
            "beast.experimental.disableBlockDiagonalNativeSelectionPath";

    private final Parameter selectionMatrixParameter;
    private final int stateDimension;
    private final OUProcessModel processModel;

    private final CanonicalBranchMessageContribution localContribution;
    private final CanonicalBranchMessageContributionUtils.Workspace contributionWorkspace;
    private final CanonicalLocalTransitionAdjoints localAdjoints;
    private final CanonicalTransitionAdjointUtils.Workspace canonicalAdjointWorkspace;
    private final double[][] dLogL_dF;
    private final double[][] dLogL_dV;
    private final double[] nativeCompressedGradientScratch;
    // Row-major flat rotation gradient for CanonicalSelectionGradientProjector.
    private final double[] nativeRotationGradientFlat;
    private final double[] stationaryMeanScratch;
    private final CanonicalBlockDiagonalGradientCache blockDiagonalGradientCache;

    public CanonicalSelectionMatrixGradientFormula(final Parameter selectionMatrixParameter,
                                                   final int stateDimension) {
        this(null, selectionMatrixParameter, stateDimension);
    }

    public CanonicalSelectionMatrixGradientFormula(final OUProcessModel processModel,
                                                   final Parameter selectionMatrixParameter,
                                                   final int stateDimension) {
        this(processModel, selectionMatrixParameter, stateDimension, null);
    }

    public CanonicalSelectionMatrixGradientFormula(final OUProcessModel processModel,
                                                   final Parameter selectionMatrixParameter,
                                                   final int stateDimension,
                                                   final CanonicalBlockDiagonalGradientCache blockDiagonalGradientCache) {
        if (selectionMatrixParameter == null) {
            throw new IllegalArgumentException("selectionMatrixParameter must not be null");
        }
        if (stateDimension < 1) {
            throw new IllegalArgumentException("stateDimension must be at least 1");
        }
        this.selectionMatrixParameter = selectionMatrixParameter;
        this.stateDimension = stateDimension;
        this.processModel = processModel;
        this.blockDiagonalGradientCache = blockDiagonalGradientCache;

        this.localContribution = new CanonicalBranchMessageContribution(stateDimension);
        this.contributionWorkspace = new CanonicalBranchMessageContributionUtils.Workspace(stateDimension);
        this.localAdjoints = new CanonicalLocalTransitionAdjoints(stateDimension);
        this.canonicalAdjointWorkspace = new CanonicalTransitionAdjointUtils.Workspace(stateDimension);
        this.dLogL_dF = new double[stateDimension][stateDimension];
        this.dLogL_dV = new double[stateDimension][stateDimension];
        final boolean hasBlock = selectionMatrixParameter instanceof AbstractBlockDiagonalTwoByTwoMatrixParameter;
        this.nativeCompressedGradientScratch = hasBlock
                ? new double[((AbstractBlockDiagonalTwoByTwoMatrixParameter) selectionMatrixParameter).getCompressedDDimension()]
                : null;
        this.nativeRotationGradientFlat = hasBlock
                ? new double[stateDimension * stateDimension]
                : null;
        this.stationaryMeanScratch = processModel != null ? new double[stateDimension] : null;
    }

    @Override
    public boolean supportsParameter(final Parameter parameter) {
        if (parameter == selectionMatrixParameter) {
            return true;
        }
        if (!(selectionMatrixParameter instanceof AbstractBlockDiagonalTwoByTwoMatrixParameter)) {
            return false;
        }

        final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter =
                (AbstractBlockDiagonalTwoByTwoMatrixParameter) selectionMatrixParameter;
        if (parameter == blockParameter.getParameter()) {
            return true;
        }
        if (parameter == blockParameter.getRotationMatrixParameter()) {
            return true;
        }
        if (blockParameter.getRotationMatrixParameter() instanceof OrthogonalMatrixProvider
                && parameter == ((OrthogonalMatrixProvider) blockParameter.getRotationMatrixParameter()).getOrthogonalParameter()) {
            return true;
        }
        if (parameter == blockParameter.getScalarBlockParameter()
                && blockParameter.getScalarBlockParameter().getDimension() > 0) {
            return true;
        }
        for (int i = 0; i < blockParameter.getTwoByTwoParameterFamilyCount(); ++i) {
            if (parameter == blockParameter.getTwoByTwoBlockParameter(i)) {
                return true;
            }
        }
        return false;
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
        if (shouldUseBlockDiagonalNativePath(parameter)) {
            if (blockDiagonalGradientCache != null
                    && blockDiagonalGradientCache.supportsNativeSelectionParameter(parameter)
                    && branchGradientCache != null) {
                return blockDiagonalGradientCache.getNativeSelectionGradient(
                        parameter, trajectory, branchGradientCache, repr, timeGrid);
            }
            return computeBlockDiagonalNativeGradient(parameter, trajectory, branchGradientCache, timeGrid);
        }

        final int d = stateDimension;
        final int T = trajectory.timeCount;
        final double[] gradientAccumulator = new double[d * d];
        if (branchGradientCache != null) {
            branchGradientCache.ensure(trajectory);
        }

        for (int t = 0; t < T - 1; ++t) {
            final CanonicalLocalTransitionAdjoints adjoints =
                    localAdjoints(t, trajectory, branchGradientCache);
            GaussianMatrixOps.copyFlatToMatrix(adjoints.dLogL_dF, dLogL_dF, d);
            repr.accumulateSelectionGradient(
                    t, t + 1, timeGrid, dLogL_dF, adjoints.dLogL_df, gradientAccumulator);

            GaussianMatrixOps.copyFlatToMatrix(adjoints.dLogL_dOmega, dLogL_dV, d);
            repr.accumulateSelectionGradientFromCovariance(t, t + 1, timeGrid, dLogL_dV, gradientAccumulator);
        }

        return gradientAccumulator;
    }

    private boolean shouldUseBlockDiagonalNativePath(final Parameter requestedParameter) {
        if (Boolean.getBoolean(DISABLE_BLOCK_DIAGONAL_NATIVE_SELECTION_PATH_PROPERTY)
                || Boolean.getBoolean(DISABLE_ORTHOGONAL_NATIVE_SELECTION_PATH_PROPERTY)) {
            return false;
        }
        if (requestedParameter == selectionMatrixParameter) {
            return false;
        }
        if (!(selectionMatrixParameter instanceof AbstractBlockDiagonalTwoByTwoMatrixParameter)) {
            return false;
        }
        return processModel != null
                && processModel.getCovarianceGradientMethod() == OUProcessModel.CovarianceGradientMethod.STATIONARY_LYAPUNOV
                && processModel.getSelectionMatrixParameterization()
                instanceof BlockDiagonalNativeCanonicalParameterization;
    }

    @Override
    public void makeDirty() {
        if (blockDiagonalGradientCache != null) {
            blockDiagonalGradientCache.makeDirty();
        }
    }

    private double[] computeBlockDiagonalNativeGradient(final Parameter requestedParameter,
                                                       final CanonicalForwardTrajectory trajectory,
                                                       final CanonicalBranchGradientCache branchGradientCache,
                                                       final TimeGrid timeGrid) {
        final BlockDiagonalNativeCanonicalParameterization blockParameterization =
                (BlockDiagonalNativeCanonicalParameterization) processModel.getSelectionMatrixParameterization();
        final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter =
                (AbstractBlockDiagonalTwoByTwoMatrixParameter) selectionMatrixParameter;
        final int T = trajectory.timeCount;
        final double[] compressedDGradient = nativeCompressedGradientScratch;
        final double[] rotationGradientFlat = nativeRotationGradientFlat;
        final double[] stationaryMean = stationaryMeanScratch;
        Arrays.fill(compressedDGradient, 0.0);
        Arrays.fill(rotationGradientFlat, 0.0);
        processModel.getInitialMean(stationaryMean);
        if (branchGradientCache != null) {
            branchGradientCache.ensure(trajectory);
        }

        for (int t = 0; t < T - 1; ++t) {
            final double dt = timeGrid.getDelta(t, t + 1);
            final CanonicalBranchMessageContribution contribution =
                    localContribution(t, trajectory, branchGradientCache);
            blockParameterization.accumulateNativeGradientFromCanonicalContributionFlat(
                    processModel.getDiffusionMatrix(),
                    stationaryMean,
                    dt,
                    contribution,
                    localAdjoints(t, trajectory, branchGradientCache),
                    compressedDGradient,
                    rotationGradientFlat);
        }

        final double[] nativeGradient = new double[blockParameter.getBlockDiagonalNParameters()];
        blockParameter.chainGradient(compressedDGradient, nativeGradient);
        return CanonicalSelectionGradientProjector.assembleBlockGradientResultFlat(
                stateDimension,
                requestedParameter,
                blockParameter,
                nativeGradient,
                rotationGradientFlat);
    }

    private CanonicalBranchMessageContribution localContribution(final int branchIndex,
                                                                 final CanonicalForwardTrajectory trajectory,
                                                                 final CanonicalBranchGradientCache branchGradientCache) {
        if (branchGradientCache != null) {
            return branchGradientCache.getContribution(branchIndex);
        }
        CanonicalBranchMessageContributionUtils.fillFromPairState(
                trajectory.branchPairStates[branchIndex],
                contributionWorkspace,
                localContribution);
        return localContribution;
    }

    private CanonicalLocalTransitionAdjoints localAdjoints(final int branchIndex,
                                                           final CanonicalForwardTrajectory trajectory,
                                                           final CanonicalBranchGradientCache branchGradientCache) {
        if (branchGradientCache != null) {
            return branchGradientCache.getAdjoints(branchIndex);
        }
        CanonicalTransitionAdjointUtils.fillFromCanonicalTransition(
                trajectory.transitions[branchIndex],
                localContribution(branchIndex, trajectory, null),
                canonicalAdjointWorkspace,
                localAdjoints);
        return localAdjoints;
    }
}
