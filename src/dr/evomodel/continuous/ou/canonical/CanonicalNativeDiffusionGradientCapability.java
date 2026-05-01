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
}
