package dr.evomodel.continuous.ou.canonical;

/**
 * Capability for native prepared-branch diffusion-gradient pullbacks.
 */
public interface CanonicalNativeDiffusionGradientCapability {

    void accumulateDiffusionGradientPreparedFlat(CanonicalPreparedBranchHandle prepared,
                                                 double[] dLogL_dV,
                                                 boolean transposeAdjoint,
                                                 double[] gradientAccumulator,
                                                 CanonicalBranchWorkspace workspace);

    default boolean supportsDelayedDiffusionGradientRotation() {
        return false;
    }

    default void accumulateDiffusionGradientPreparedDBasisFlat(CanonicalPreparedBranchHandle prepared,
                                                               double[] dLogL_dV,
                                                               boolean transposeAdjoint,
                                                               double[] dBasisGradientAccumulator,
                                                               CanonicalBranchWorkspace workspace) {
        throw new UnsupportedOperationException(
                "Delayed diffusion-gradient rotation is not supported by this canonical OU backend.");
    }

    default void finishDiffusionGradientFromDBasisFlat(double[] dBasisGradientAccumulator,
                                                       double[] gradientAccumulator,
                                                       CanonicalBranchWorkspace workspace) {
        throw new UnsupportedOperationException(
                "Delayed diffusion-gradient rotation is not supported by this canonical OU backend.");
    }
}
