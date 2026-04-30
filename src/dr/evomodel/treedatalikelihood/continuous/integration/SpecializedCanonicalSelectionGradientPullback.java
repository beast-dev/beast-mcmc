package dr.evomodel.treedatalikelihood.continuous.integration;

import dr.evomodel.continuous.ou.CanonicalBranchWorkspace;
import dr.evomodel.continuous.ou.CanonicalPreparedBranchHandle;
import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.continuous.ou.SpecializedCanonicalSelectionParameterization;

import java.util.Arrays;

final class SpecializedCanonicalSelectionGradientPullback implements CanonicalSelectionGradientPullback {

    private final SpecializedCanonicalSelectionParameterization selection;
    private final int dimension;
    private final int compressedGradientDim;
    private final int nativeGradientScratchDim;

    SpecializedCanonicalSelectionGradientPullback(
            final SpecializedCanonicalSelectionParameterization selection,
            final int dimension,
            final double[] gradA,
            final BranchGradientWorkspace workspace) {
        this.selection = selection;
        this.dimension = dimension;
        if (gradA.length != selection.getSelectionGradientDimension()) {
            throw new IllegalArgumentException(
                    "Specialized selection gradient expects native parameter length "
                            + selection.getSelectionGradientDimension() + ", found " + gradA.length);
        }
        this.compressedGradientDim = selection.getCompressedSelectionGradientDimension();
        this.nativeGradientScratchDim = selection.getNativeSelectionGradientScratchDimension();
        final GradientPullbackWorkspace gradient = workspace.gradient;
        if (compressedGradientDim > gradient.orthogonalCompressedGradientScratch.length
                || nativeGradientScratchDim > gradient.orthogonalNativeGradientScratch.length) {
            throw new IllegalStateException(
                    "Specialized selection scratch is too small for native gradient dimensions "
                            + compressedGradientDim + " and " + nativeGradientScratchDim + ".");
        }
    }

    @Override
    public void initialize(final BranchGradientWorkspace workspace,
                           final double[] gradA,
                           final double[] gradMu) {
        final GradientPullbackWorkspace gradient = workspace.gradient;
        Arrays.fill(gradient.orthogonalCompressedGradientScratch, 0, compressedGradientDim, 0.0);
        Arrays.fill(gradient.orthogonalNativeGradientScratch, 0, nativeGradientScratchDim, 0.0);
        Arrays.fill(gradient.orthogonalRotationGradientFlatScratch, 0.0);
        workspace.ensureSpecializedBranchWorkspace(selection);
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
                compressedGradientDim);
    }

    @Override
    public void prepareWorkspace(final BranchGradientWorkspace workspace) {
        workspace.ensureSpecializedBranchWorkspace(selection);
    }

    @Override
    public void accumulateForBranch(final OUProcessModel processModel,
                                    final BranchGradientInputs inputs,
                                    final int activeIndex,
                                    final BranchGradientWorkspace workspace,
                                    final double[] gradA,
                                    final double[] gradQ,
                                    final double[] gradMu) {
        final CanonicalPreparedBranchHandle prepared = inputs.getPreparedBranchHandle(activeIndex);
        final CanonicalBranchWorkspace specializedWorkspace =
                workspace.ensureSpecializedBranchWorkspace(selection);
        final GradientPullbackWorkspace gradient = workspace.gradient;

        selection.accumulateNativeGradientFromAdjointsPreparedFlat(
                prepared,
                processModel.getDiffusionMatrix(),
                workspace.adjoints,
                specializedWorkspace,
                gradient.orthogonalCompressedGradientScratch,
                gradient.orthogonalRotationGradientFlatScratch);

        selection.accumulateDiffusionGradientPreparedFlat(
                prepared,
                workspace.adjoints.dLogL_dOmega,
                true,
                gradQ,
                specializedWorkspace);

        selection.accumulateMeanGradientPrepared(
                prepared,
                workspace.adjoints.dLogL_df,
                gradMu,
                specializedWorkspace);
    }

    @Override
    public void reduceWorker(final BranchGradientWorkspace worker,
                             final BranchGradientWorkspace reductionWorkspace,
                             final double[] gradA) {
        accumulateVectorInPlace(
                reductionWorkspace.gradient.orthogonalCompressedGradientScratch,
                worker.gradient.orthogonalCompressedGradientScratch,
                compressedGradientDim);
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
        selection.finishNativeSelectionGradient(
                gradient.orthogonalCompressedGradientScratch,
                gradient.orthogonalNativeGradientScratch,
                gradient.orthogonalRotationGradientFlatScratch,
                gradA);
    }

    private static void accumulateVectorInPlace(final double[] target,
                                                final double[] source,
                                                final int length) {
        for (int i = 0; i < length; ++i) {
            target[i] += source[i];
        }
    }
}
