package dr.evomodel.continuous.ou.canonical;

/**
 * Capability for native prepared-branch stationary-mean pullbacks.
 */
public interface CanonicalNativeMeanGradientCapability {

    void accumulateMeanGradientPrepared(CanonicalPreparedBranchHandle prepared,
                                        double[] dLogL_df,
                                        double[] gradientAccumulator,
                                        CanonicalBranchWorkspace workspace);

    default boolean supportsDelayedMeanGradientRotation() {
        return false;
    }

    default void accumulateMeanGradientPreparedDBasisFlat(CanonicalPreparedBranchHandle prepared,
                                                          double[] dLogL_df,
                                                          double[] dBasisGradientAccumulator,
                                                          CanonicalBranchWorkspace workspace) {
        throw new UnsupportedOperationException(
                "Delayed stationary-mean gradient rotation is not supported by this canonical OU backend.");
    }

    default void finishMeanGradientFromDBasisFlat(double[] dBasisGradientAccumulator,
                                                  double[] gradientAccumulator,
                                                  CanonicalBranchWorkspace workspace) {
        throw new UnsupportedOperationException(
                "Delayed stationary-mean gradient rotation is not supported by this canonical OU backend.");
    }
}
