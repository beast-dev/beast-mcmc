package dr.evomodel.treedatalikelihood.continuous.canonical.gradient;

import dr.evomodel.continuous.ou.canonical.CanonicalBranchWorkspace;
import dr.evomodel.continuous.ou.canonical.CanonicalPreparedBranchGradientAccumulator;
import dr.evomodel.continuous.ou.canonical.CanonicalPreparedBranchHandle;
import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.continuous.ou.canonical.SpecializedCanonicalSelectionParameterization;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalLocalTransitionAdjoints;
import dr.evomodel.treedatalikelihood.continuous.canonical.workspace.BranchGradientWorkspace;
import dr.evomodel.treedatalikelihood.continuous.canonical.workspace.SpecializedGradientWorkspace;

final class SpecializedCanonicalSelectionGradientPullback implements CanonicalSelectionGradientPullback {

    private final SpecializedCanonicalSelectionParameterization selection;
    private final int dimension;
    private final CanonicalPreparedBranchGradientAccumulator accumulator;

    SpecializedCanonicalSelectionGradientPullback(
            final SpecializedCanonicalSelectionParameterization selection,
            final int dimension,
            final double[] gradA,
            final BranchGradientWorkspace workspace) {
        this.selection = selection;
        this.dimension = dimension;
        this.accumulator = new CanonicalPreparedBranchGradientAccumulator(selection, dimension);
        if (gradA.length != selection.getSelectionGradientDimension()) {
            throw new IllegalArgumentException(
                    "Specialized selection gradient expects native parameter length "
                            + selection.getSelectionGradientDimension() + ", found " + gradA.length);
        }
        final SpecializedGradientWorkspace gradient = workspace.specializedGradient();
        accumulator.checkBuffers(
                gradient.compressedGradientScratch,
                gradient.nativeGradientScratch,
                gradient.rotationGradientFlatScratch,
                gradient.diffusionGradientDBasisScratch,
                gradient.meanGradientDBasisScratch);
    }

    @Override
    public void initialize(final BranchGradientWorkspace workspace,
                           final double[] gradA,
                           final double[] gradMu) {
        final SpecializedGradientWorkspace gradient = workspace.specializedGradient();
        accumulator.clearBuffers(
                gradient.compressedGradientScratch,
                gradient.nativeGradientScratch,
                gradient.rotationGradientFlatScratch,
                gradient.diffusionGradientDBasisScratch,
                gradient.meanGradientDBasisScratch);
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
                accumulator.getCompressedSelectionGradientDimension());
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

        accumulator.accumulateSelectionDiffusionAndMean(
                prepared,
                processModel.getDiffusionMatrix(),
                localAdjoints,
                specializedWorkspace,
                gradient.compressedGradientScratch,
                gradient.rotationGradientFlatScratch,
                gradQ,
                gradient.diffusionGradientDBasisScratch,
                gradMu,
                gradient.meanGradientDBasisScratch);
    }

    @Override
    public void reduceWorker(final BranchGradientWorkspace worker,
                             final BranchGradientWorkspace reductionWorkspace,
                             final double[] gradA) {
        accumulateVectorInPlace(
                reductionWorkspace.specializedGradient().compressedGradientScratch,
                worker.specializedGradient().compressedGradientScratch,
                accumulator.getCompressedSelectionGradientDimension());
        accumulateVectorInPlace(
                reductionWorkspace.specializedGradient().rotationGradientFlatScratch,
                worker.specializedGradient().rotationGradientFlatScratch,
                reductionWorkspace.specializedGradient().rotationGradientFlatScratch.length);
        if (accumulator.delaysDiffusionGradientRotation()) {
            accumulateVectorInPlace(
                    reductionWorkspace.specializedGradient().diffusionGradientDBasisScratch,
                    worker.specializedGradient().diffusionGradientDBasisScratch,
                    reductionWorkspace.specializedGradient().diffusionGradientDBasisScratch.length);
        }
        if (accumulator.delaysMeanGradientRotation()) {
            accumulateVectorInPlace(
                    reductionWorkspace.specializedGradient().meanGradientDBasisScratch,
                    worker.specializedGradient().meanGradientDBasisScratch,
                    reductionWorkspace.specializedGradient().meanGradientDBasisScratch.length);
        }
    }

    @Override
    public void finish(final BranchGradientInputs inputs,
                       final BranchGradientWorkspace workspace,
                       final double[] gradA,
                       final double[] gradQ,
                       final double[] gradMu) {
        final SpecializedGradientWorkspace gradient = workspace.specializedGradient();
        final CanonicalBranchWorkspace specializedWorkspace =
                workspace.ensureSpecializedBranchWorkspace(selection);
        accumulator.finishSelection(
                gradient.compressedGradientScratch,
                gradient.nativeGradientScratch,
                gradient.rotationGradientFlatScratch,
                gradA);
        accumulator.finishDelayedDiffusion(
                gradient.diffusionGradientDBasisScratch,
                gradQ,
                specializedWorkspace);
        accumulator.finishDelayedMean(
                gradient.meanGradientDBasisScratch,
                gradMu,
                specializedWorkspace);
    }

    private static void accumulateVectorInPlace(final double[] target,
                                                final double[] source,
                                                final int length) {
        for (int i = 0; i < length; ++i) {
            target[i] += source[i];
        }
    }
}
