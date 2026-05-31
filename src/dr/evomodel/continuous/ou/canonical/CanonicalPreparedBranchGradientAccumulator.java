package dr.evomodel.continuous.ou.canonical;

import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalLocalTransitionAdjoints;
import dr.inference.model.MatrixParameterInterface;

import java.util.Arrays;

/**
 * Shared prepared-branch pullback for canonical OU gradients.
 *
 * <p>The tree and time-series canonical paths both reduce their local factors
 * to {@link CanonicalLocalTransitionAdjoints}. This class owns the common
 * native prepared-branch accumulation from those adjoints to selection,
 * diffusion, and stationary-mean gradient buffers.</p>
 */
public final class CanonicalPreparedBranchGradientAccumulator {

    private final SpecializedCanonicalSelectionParameterization selection;
    private final int dimension;
    private final int compressedSelectionGradientDimension;
    private final int nativeSelectionGradientScratchDimension;
    private final boolean delayDiffusionGradientRotation;
    private final boolean delayMeanGradientRotation;

    public CanonicalPreparedBranchGradientAccumulator(
            final SpecializedCanonicalSelectionParameterization selection,
            final int dimension) {
        if (selection == null) {
            throw new IllegalArgumentException("selection must not be null");
        }
        if (dimension < 1) {
            throw new IllegalArgumentException("dimension must be positive");
        }
        this.selection = selection;
        this.dimension = dimension;
        this.compressedSelectionGradientDimension =
                selection.getCompressedSelectionGradientDimension();
        this.nativeSelectionGradientScratchDimension =
                selection.getNativeSelectionGradientScratchDimension();
        this.delayDiffusionGradientRotation =
                selection.supportsDelayedDiffusionGradientRotation();
        this.delayMeanGradientRotation =
                selection.supportsDelayedMeanGradientRotation();
    }

    public int getCompressedSelectionGradientDimension() {
        return compressedSelectionGradientDimension;
    }

    public int getNativeSelectionGradientScratchDimension() {
        return nativeSelectionGradientScratchDimension;
    }

    public boolean delaysDiffusionGradientRotation() {
        return delayDiffusionGradientRotation;
    }

    public boolean delaysMeanGradientRotation() {
        return delayMeanGradientRotation;
    }

    public void checkBuffers(final double[] compressedSelectionGradient,
                             final double[] nativeSelectionGradientScratch,
                             final double[] rotationGradientFlat,
                             final double[] diffusionGradientDBasis,
                             final double[] meanGradientDBasis) {
        requireLength(
                compressedSelectionGradient,
                compressedSelectionGradientDimension,
                "compressedSelectionGradient");
        requireLength(
                nativeSelectionGradientScratch,
                nativeSelectionGradientScratchDimension,
                "nativeSelectionGradientScratch");
        requireLength(rotationGradientFlat, dimension * dimension, "rotationGradientFlat");
        requireLength(diffusionGradientDBasis, dimension * dimension, "diffusionGradientDBasis");
        requireLength(meanGradientDBasis, dimension, "meanGradientDBasis");
    }

    public void clearBuffers(final double[] compressedSelectionGradient,
                             final double[] nativeSelectionGradientScratch,
                             final double[] rotationGradientFlat,
                             final double[] diffusionGradientDBasis,
                             final double[] meanGradientDBasis) {
        checkBuffers(
                compressedSelectionGradient,
                nativeSelectionGradientScratch,
                rotationGradientFlat,
                diffusionGradientDBasis,
                meanGradientDBasis);
        Arrays.fill(compressedSelectionGradient, 0, compressedSelectionGradientDimension, 0.0);
        Arrays.fill(nativeSelectionGradientScratch, 0, nativeSelectionGradientScratchDimension, 0.0);
        Arrays.fill(rotationGradientFlat, 0.0);
        Arrays.fill(diffusionGradientDBasis, 0.0);
        Arrays.fill(meanGradientDBasis, 0.0);
    }

    public void accumulateSelectionAndMean(final CanonicalPreparedBranchHandle prepared,
                                           final MatrixParameterInterface diffusionMatrix,
                                           final CanonicalLocalTransitionAdjoints localAdjoints,
                                           final CanonicalBranchWorkspace workspace,
                                           final double[] compressedSelectionGradient,
                                           final double[] rotationGradientFlat,
                                           final double[] meanGradient,
                                           final double[] meanGradientDBasis) {
        selection.accumulateNativeGradientFromAdjointsPreparedFlat(
                prepared,
                diffusionMatrix,
                localAdjoints,
                workspace,
                compressedSelectionGradient,
                rotationGradientFlat);
        accumulateMean(prepared, localAdjoints, workspace, meanGradient, meanGradientDBasis);
    }

    public void accumulateDiffusion(final CanonicalPreparedBranchHandle prepared,
                                    final CanonicalLocalTransitionAdjoints localAdjoints,
                                    final CanonicalBranchWorkspace workspace,
                                    final boolean transposeDiffusionAdjoint,
                                    final double[] diffusionGradient,
                                    final double[] diffusionGradientDBasis) {
        accumulateDiffusion(
                prepared,
                localAdjoints,
                workspace,
                transposeDiffusionAdjoint,
                true,
                diffusionGradient,
                diffusionGradientDBasis);
    }

    public void accumulateDiffusion(final CanonicalPreparedBranchHandle prepared,
                                    final CanonicalLocalTransitionAdjoints localAdjoints,
                                    final CanonicalBranchWorkspace workspace,
                                    final boolean transposeDiffusionAdjoint,
                                    final boolean allowDelayedDiffusionRotation,
                                    final double[] diffusionGradient,
                                    final double[] diffusionGradientDBasis) {
        if (allowDelayedDiffusionRotation && delayDiffusionGradientRotation) {
            selection.accumulateDiffusionGradientPreparedDBasisFlat(
                    prepared,
                    localAdjoints.dLogL_dOmega,
                    transposeDiffusionAdjoint,
                    diffusionGradientDBasis,
                    workspace);
        } else {
            selection.accumulateDiffusionGradientPreparedFlat(
                    prepared,
                    localAdjoints.dLogL_dOmega,
                    transposeDiffusionAdjoint,
                    diffusionGradient,
                    workspace);
        }
    }

    public void accumulateSelectionDiffusionAndMean(
            final CanonicalPreparedBranchHandle prepared,
            final MatrixParameterInterface diffusionMatrix,
            final CanonicalLocalTransitionAdjoints localAdjoints,
            final CanonicalBranchWorkspace workspace,
            final double[] compressedSelectionGradient,
            final double[] rotationGradientFlat,
            final double[] diffusionGradient,
            final double[] diffusionGradientDBasis,
            final double[] meanGradient,
            final double[] meanGradientDBasis) {
        selection.accumulateNativeSelectionAndDiffusionGradientFromAdjointsPreparedFlat(
                prepared,
                diffusionMatrix,
                localAdjoints,
                workspace,
                compressedSelectionGradient,
                rotationGradientFlat,
                delayDiffusionGradientRotation,
                delayDiffusionGradientRotation ? diffusionGradientDBasis : diffusionGradient);
        accumulateMean(prepared, localAdjoints, workspace, meanGradient, meanGradientDBasis);
    }

    public void finishSelection(final double[] compressedSelectionGradient,
                                final double[] nativeSelectionGradientScratch,
                                final double[] rotationGradientFlat,
                                final double[] selectionGradient) {
        selection.finishNativeSelectionGradient(
                compressedSelectionGradient,
                nativeSelectionGradientScratch,
                rotationGradientFlat,
                selectionGradient);
    }

    public void finishDelayedDiffusion(final double[] diffusionGradientDBasis,
                                       final double[] diffusionGradient,
                                       final CanonicalBranchWorkspace workspace) {
        if (!delayDiffusionGradientRotation) {
            return;
        }
        selection.finishDiffusionGradientFromDBasisFlat(
                diffusionGradientDBasis,
                diffusionGradient,
                workspace);
    }

    public void finishDelayedMean(final double[] meanGradientDBasis,
                                  final double[] meanGradient,
                                  final CanonicalBranchWorkspace workspace) {
        if (!delayMeanGradientRotation) {
            return;
        }
        selection.finishMeanGradientFromDBasisFlat(
                meanGradientDBasis,
                meanGradient,
                workspace);
    }

    private void accumulateMean(final CanonicalPreparedBranchHandle prepared,
                                final CanonicalLocalTransitionAdjoints localAdjoints,
                                final CanonicalBranchWorkspace workspace,
                                final double[] meanGradient,
                                final double[] meanGradientDBasis) {
        if (delayMeanGradientRotation) {
            selection.accumulateMeanGradientPreparedDBasisFlat(
                    prepared,
                    localAdjoints.dLogL_df,
                    meanGradientDBasis,
                    workspace);
        } else {
            selection.accumulateMeanGradientPrepared(
                    prepared,
                    localAdjoints.dLogL_df,
                    meanGradient,
                    workspace);
        }
    }

    private static void requireLength(final double[] array,
                                      final int expectedLength,
                                      final String name) {
        if (array == null || array.length < expectedLength) {
            throw new IllegalArgumentException(
                    name + " length must be at least " + expectedLength);
        }
    }
}
