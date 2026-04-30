package dr.evomodel.continuous.ou;

/**
 * Capability for native prepared-branch stationary-mean pullbacks.
 */
public interface CanonicalNativeMeanGradientCapability {

    void accumulateMeanGradientPrepared(CanonicalPreparedBranchHandle prepared,
                                        double[] dLogL_df,
                                        double[] gradientAccumulator,
                                        CanonicalBranchWorkspace workspace);
}
