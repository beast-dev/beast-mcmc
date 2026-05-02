package dr.evomodel.continuous.ou.canonical;

import dr.evomodel.continuous.ou.SelectionMatrixParameterization;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalLocalTransitionAdjoints;
import dr.inference.model.MatrixParameterInterface;

/**
 * Optional canonical OU backend for selection parametrizations with a faster
 * native transition, covariance, or gradient path.
 */
public interface SpecializedCanonicalSelectionParameterization extends SelectionMatrixParameterization,
        CanonicalGradientPackingCapability,
        CanonicalPreparedTransitionCapability,
        CanonicalNativeSelectionGradientCapability,
        CanonicalNativeBranchGradientCapability,
        CanonicalNativeMeanGradientCapability,
        CanonicalNativeDiffusionGradientCapability {

    default void accumulateNativeSelectionAndDiffusionGradientFromAdjointsPreparedFlat(
            final CanonicalPreparedBranchHandle prepared,
            final MatrixParameterInterface diffusionMatrix,
            final CanonicalLocalTransitionAdjoints localAdjoints,
            final CanonicalBranchWorkspace workspace,
            final double[] compressedDAccumulator,
            final double[] rotationAccumulator,
            final boolean delayDiffusionGradientRotation,
            final double[] diffusionGradientAccumulator) {
        accumulateNativeGradientFromAdjointsPreparedFlat(
                prepared,
                diffusionMatrix,
                localAdjoints,
                workspace,
                compressedDAccumulator,
                rotationAccumulator);

        if (delayDiffusionGradientRotation) {
            accumulateDiffusionGradientPreparedDBasisFlat(
                    prepared,
                    localAdjoints.dLogL_dOmega,
                    true,
                    diffusionGradientAccumulator,
                    workspace);
        } else {
            accumulateDiffusionGradientPreparedFlat(
                    prepared,
                    localAdjoints.dLogL_dOmega,
                    true,
                    diffusionGradientAccumulator,
                    workspace);
        }
    }
}
