package dr.evomodel.treedatalikelihood.continuous.canonical.gradient;

import dr.evomodel.continuous.ou.canonical.CanonicalBranchWorkspace;
import dr.evomodel.continuous.ou.canonical.CanonicalPreparedBranchHandle;
import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.continuous.ou.canonical.SpecializedCanonicalSelectionParameterization;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalLocalTransitionAdjoints;
import dr.evomodel.treedatalikelihood.continuous.canonical.workspace.BranchGradientWorkspace;
import dr.evomodel.treedatalikelihood.continuous.canonical.workspace.SpecializedGradientWorkspace;

import java.util.Arrays;

final class SpecializedCanonicalSelectionGradientPullback implements CanonicalSelectionGradientPullback {

    private final SpecializedCanonicalSelectionParameterization selection;
    private final int dimension;
    private final int compressedGradientDim;
    private final int nativeGradientScratchDim;
    private final boolean delayDiffusionGradientRotation;

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
        this.delayDiffusionGradientRotation = selection.supportsDelayedDiffusionGradientRotation();
        final SpecializedGradientWorkspace gradient = workspace.specializedGradient();
        if (compressedGradientDim > gradient.compressedGradientScratch.length
                || nativeGradientScratchDim > gradient.nativeGradientScratch.length) {
            throw new IllegalStateException(
                    "Specialized selection scratch is too small for native gradient dimensions "
                            + compressedGradientDim + " and " + nativeGradientScratchDim + ".");
        }
    }

    @Override
    public void initialize(final BranchGradientWorkspace workspace,
                           final double[] gradA,
                           final double[] gradMu) {
        final SpecializedGradientWorkspace gradient = workspace.specializedGradient();
        Arrays.fill(gradient.compressedGradientScratch, 0, compressedGradientDim, 0.0);
        Arrays.fill(gradient.nativeGradientScratch, 0, nativeGradientScratchDim, 0.0);
        Arrays.fill(gradient.rotationGradientFlatScratch, 0.0);
        Arrays.fill(gradient.diffusionGradientDBasisScratch, 0.0);
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
                                    final CanonicalLocalTransitionAdjoints localAdjoints,
                                    final BranchGradientWorkspace workspace,
                                    final double[] gradA,
                                    final double[] gradQ,
                                    final double[] gradMu) {
        final CanonicalPreparedBranchHandle prepared = inputs.getPreparedBranchHandle(activeIndex);
        final CanonicalBranchWorkspace specializedWorkspace =
                workspace.ensureSpecializedBranchWorkspace(selection);
        final SpecializedGradientWorkspace gradient = workspace.specializedGradient();

        selection.accumulateNativeGradientFromAdjointsPreparedFlat(
                prepared,
                processModel.getDiffusionMatrix(),
                localAdjoints,
                specializedWorkspace,
                gradient.compressedGradientScratch,
                gradient.rotationGradientFlatScratch);

        if (delayDiffusionGradientRotation) {
            selection.accumulateDiffusionGradientPreparedDBasisFlat(
                    prepared,
                    localAdjoints.dLogL_dOmega,
                    true,
                    gradient.diffusionGradientDBasisScratch,
                    specializedWorkspace);
        } else {
            selection.accumulateDiffusionGradientPreparedFlat(
                    prepared,
                    localAdjoints.dLogL_dOmega,
                    true,
                    gradQ,
                    specializedWorkspace);
        }

        selection.accumulateMeanGradientPrepared(
                prepared,
                localAdjoints.dLogL_df,
                gradMu,
                specializedWorkspace);
    }

    @Override
    public void reduceWorker(final BranchGradientWorkspace worker,
                             final BranchGradientWorkspace reductionWorkspace,
                             final double[] gradA) {
        accumulateVectorInPlace(
                reductionWorkspace.specializedGradient().compressedGradientScratch,
                worker.specializedGradient().compressedGradientScratch,
                compressedGradientDim);
        accumulateVectorInPlace(
                reductionWorkspace.specializedGradient().rotationGradientFlatScratch,
                worker.specializedGradient().rotationGradientFlatScratch,
                reductionWorkspace.specializedGradient().rotationGradientFlatScratch.length);
        if (delayDiffusionGradientRotation) {
            accumulateVectorInPlace(
                    reductionWorkspace.specializedGradient().diffusionGradientDBasisScratch,
                    worker.specializedGradient().diffusionGradientDBasisScratch,
                    reductionWorkspace.specializedGradient().diffusionGradientDBasisScratch.length);
        }
    }

    @Override
    public void finish(final BranchGradientInputs inputs,
                       final BranchGradientWorkspace workspace,
                       final double[] gradA,
                       final double[] gradQ) {
        final SpecializedGradientWorkspace gradient = workspace.specializedGradient();
        selection.finishNativeSelectionGradient(
                gradient.compressedGradientScratch,
                gradient.nativeGradientScratch,
                gradient.rotationGradientFlatScratch,
                gradA);
        if (delayDiffusionGradientRotation) {
            selection.finishDiffusionGradientFromDBasisFlat(
                    gradient.diffusionGradientDBasisScratch,
                    gradQ,
                    workspace.ensureSpecializedBranchWorkspace(selection));
        }
    }

    private static void accumulateVectorInPlace(final double[] target,
                                                final double[] source,
                                                final int length) {
        for (int i = 0; i < length; ++i) {
            target[i] += source[i];
        }
    }
}
