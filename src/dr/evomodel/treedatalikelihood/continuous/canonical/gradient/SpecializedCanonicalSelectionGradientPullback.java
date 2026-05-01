package dr.evomodel.treedatalikelihood.continuous.canonical.gradient;

import dr.evomodel.continuous.ou.canonical.CanonicalBranchWorkspace;
import dr.evomodel.continuous.ou.canonical.CanonicalPreparedBranchHandle;
import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.continuous.ou.canonical.SpecializedCanonicalSelectionParameterization;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalLocalTransitionAdjoints;
import dr.evomodel.treedatalikelihood.continuous.canonical.workspace.BranchGradientWorkspace;
import dr.evomodel.treedatalikelihood.continuous.canonical.workspace.OrthogonalBlockGradientWorkspace;

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
        final OrthogonalBlockGradientWorkspace gradient = workspace.orthogonalGradient();
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
        final OrthogonalBlockGradientWorkspace gradient = workspace.orthogonalGradient();
        Arrays.fill(gradient.compressedGradientScratch, 0, compressedGradientDim, 0.0);
        Arrays.fill(gradient.nativeGradientScratch, 0, nativeGradientScratchDim, 0.0);
        Arrays.fill(gradient.rotationGradientFlatScratch, 0.0);
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
        final OrthogonalBlockGradientWorkspace gradient = workspace.orthogonalGradient();

        selection.accumulateNativeGradientFromAdjointsPreparedFlat(
                prepared,
                processModel.getDiffusionMatrix(),
                localAdjoints,
                specializedWorkspace,
                gradient.compressedGradientScratch,
                gradient.rotationGradientFlatScratch);

        selection.accumulateDiffusionGradientPreparedFlat(
                prepared,
                localAdjoints.dLogL_dOmega,
                true,
                gradQ,
                specializedWorkspace);

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
                reductionWorkspace.orthogonalGradient().compressedGradientScratch,
                worker.orthogonalGradient().compressedGradientScratch,
                compressedGradientDim);
        accumulateVectorInPlace(
                reductionWorkspace.orthogonalGradient().rotationGradientFlatScratch,
                worker.orthogonalGradient().rotationGradientFlatScratch,
                reductionWorkspace.orthogonalGradient().rotationGradientFlatScratch.length);
    }

    @Override
    public void finish(final BranchGradientInputs inputs,
                       final BranchGradientWorkspace workspace,
                       final double[] gradA) {
        final OrthogonalBlockGradientWorkspace gradient = workspace.orthogonalGradient();
        selection.finishNativeSelectionGradient(
                gradient.compressedGradientScratch,
                gradient.nativeGradientScratch,
                gradient.rotationGradientFlatScratch,
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
