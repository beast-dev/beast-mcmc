package dr.evomodel.continuous.ou;

import dr.evomodel.treedatalikelihood.continuous.gaussian.message.CanonicalLocalTransitionAdjoints;
import dr.inference.model.MatrixParameterInterface;

/**
 * Capability for native prepared-branch selection-gradient pullbacks.
 */
public interface CanonicalNativeSelectionGradientCapability {

    void accumulateNativeGradientFromAdjointsPreparedFlat(
            CanonicalPreparedBranchHandle prepared,
            MatrixParameterInterface diffusionMatrix,
            CanonicalLocalTransitionAdjoints localAdjoints,
            CanonicalBranchWorkspace workspace,
            double[] compressedDAccumulator,
            double[] rotationAccumulator);
}
