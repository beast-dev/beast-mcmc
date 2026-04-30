package dr.evomodel.treedatalikelihood.continuous.integration;

import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.continuous.ou.orthogonalblockdiagonal.OrthogonalBlockBranchGradientWorkspace;
import dr.evomodel.continuous.ou.orthogonalblockdiagonal.OrthogonalBlockCanonicalParameterization;
import dr.evomodel.continuous.ou.orthogonalblockdiagonal.OrthogonalBlockPreparedBranchBasis;
import dr.inference.model.AbstractBlockDiagonalTwoByTwoMatrixParameter;
import dr.inference.model.OrthogonalMatrixProvider;

import java.util.Arrays;

final class OrthogonalBlockCanonicalSelectionGradientPullback implements CanonicalSelectionGradientPullback {

    private final OrthogonalBlockCanonicalParameterization orthogonalSelection;
    private final int dimension;
    private final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter;
    private final OrthogonalMatrixProvider orthogonalRotation;
    private final int compressedBlockDim;
    private final int nativeBlockDim;

    OrthogonalBlockCanonicalSelectionGradientPullback(
            final OrthogonalBlockCanonicalParameterization orthogonalSelection,
            final int dimension,
            final double[] gradA,
            final BranchGradientWorkspace workspace) {
        this.orthogonalSelection = orthogonalSelection;
        this.dimension = dimension;
        this.blockParameter =
                (AbstractBlockDiagonalTwoByTwoMatrixParameter) orthogonalSelection.getMatrixParameter();
        if (!(blockParameter.getRotationMatrixParameter() instanceof OrthogonalMatrixProvider)) {
            throw new IllegalStateException(
                    "Orthogonal block native gradient requires an OrthogonalMatrixProvider rotation parameter.");
        }
        this.orthogonalRotation =
                (OrthogonalMatrixProvider) blockParameter.getRotationMatrixParameter();
        this.nativeBlockDim = blockParameter.getBlockDiagonalNParameters();
        final int angleDim = orthogonalRotation.getOrthogonalParameter().getDimension();
        final int nativeDim = nativeBlockDim + angleDim;
        if (gradA.length != nativeDim) {
            throw new IllegalArgumentException(
                    "Orthogonal block selection gradient expects native parameter length "
                            + nativeDim + ", found " + gradA.length);
        }

        this.compressedBlockDim = blockParameter.getCompressedDDimension();
        final GradientPullbackWorkspace gradient = workspace.gradient;
        if (compressedBlockDim > gradient.orthogonalCompressedGradientScratch.length
                || nativeBlockDim > gradient.orthogonalNativeGradientScratch.length) {
            throw new IllegalStateException(
                    "Orthogonal block scratch is too small for native gradient dimensions "
                            + compressedBlockDim + " and " + nativeBlockDim + ".");
        }
    }

    @Override
    public void initialize(final BranchGradientWorkspace workspace,
                           final double[] gradA,
                           final double[] gradMu) {
        final GradientPullbackWorkspace gradient = workspace.gradient;
        Arrays.fill(gradient.orthogonalCompressedGradientScratch, 0, compressedBlockDim, 0.0);
        Arrays.fill(gradient.orthogonalNativeGradientScratch, 0, nativeBlockDim, 0.0);
        Arrays.fill(gradient.orthogonalRotationGradientFlatScratch, 0.0);
        workspace.ensureOrthogonalBranchWorkspace(orthogonalSelection);
    }

    @Override
    public void clearWorkerBuffers(final BranchGradientWorkspace workspace,
                                   final int gradALength,
                                   final int gradMuLength) {
        workspace.clearLocalGradientBuffers(
                gradALength,
                gradMuLength,
                dimension,
                true,
                compressedBlockDim);
    }

    @Override
    public void prepareWorkspace(final BranchGradientWorkspace workspace) {
        workspace.ensureOrthogonalBranchWorkspace(orthogonalSelection);
    }

    @Override
    public void accumulateForBranch(final OUProcessModel processModel,
                                    final BranchGradientInputs inputs,
                                    final int activeIndex,
                                    final BranchGradientWorkspace workspace,
                                    final double[] gradA,
                                    final double[] gradQ,
                                    final double[] gradMu) {
        final OrthogonalBlockPreparedBranchBasis preparedBasis =
                inputs.getOrthogonalPreparedBasis(activeIndex);
        final OrthogonalBlockBranchGradientWorkspace orthogonalWorkspace =
                workspace.ensureOrthogonalBranchWorkspace(orthogonalSelection);
        final GradientPullbackWorkspace gradient = workspace.gradient;

        orthogonalSelection.accumulateNativeGradientFromAdjointsPreparedFlat(
                preparedBasis,
                processModel.getDiffusionMatrix(),
                workspace.adjoints,
                orthogonalWorkspace,
                gradient.orthogonalCompressedGradientScratch,
                gradient.orthogonalRotationGradientFlatScratch);

        orthogonalSelection.accumulateDiffusionGradientPreparedFlat(
                preparedBasis,
                workspace.adjoints.dLogL_dOmega,
                true,
                gradQ,
                orthogonalWorkspace);

        orthogonalSelection.accumulateMeanGradientPrepared(
                preparedBasis,
                workspace.adjoints.dLogL_df,
                gradMu,
                orthogonalWorkspace);
    }

    @Override
    public void reduceWorker(final BranchGradientWorkspace worker,
                             final BranchGradientWorkspace reductionWorkspace,
                             final double[] gradA) {
        accumulateVectorInPlace(
                reductionWorkspace.gradient.orthogonalCompressedGradientScratch,
                worker.gradient.orthogonalCompressedGradientScratch,
                compressedBlockDim);
        accumulateVectorInPlace(
                reductionWorkspace.gradient.orthogonalRotationGradientFlatScratch,
                worker.gradient.orthogonalRotationGradientFlatScratch,
                reductionWorkspace.gradient.orthogonalRotationGradientFlatScratch.length);
    }

    @Override
    public void finish(final BranchGradientInputs inputs,
                       final BranchGradientWorkspace workspace,
                       final double[] gradA) {
        final GradientPullbackWorkspace gradient = workspace.gradient;
        blockParameter.chainGradient(
                gradient.orthogonalCompressedGradientScratch,
                gradient.orthogonalNativeGradientScratch);
        final double[] angleGradient =
                orthogonalRotation.pullBackGradientFlat(
                        gradient.orthogonalRotationGradientFlatScratch,
                        dimension);
        System.arraycopy(gradient.orthogonalNativeGradientScratch, 0, gradA, 0, nativeBlockDim);
        System.arraycopy(angleGradient, 0, gradA, nativeBlockDim, angleGradient.length);
    }

    private static void accumulateVectorInPlace(final double[] target,
                                                final double[] source,
                                                final int length) {
        for (int i = 0; i < length; ++i) {
            target[i] += source[i];
        }
    }
}
