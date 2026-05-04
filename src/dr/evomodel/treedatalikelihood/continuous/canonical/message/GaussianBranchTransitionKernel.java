package dr.evomodel.treedatalikelihood.continuous.canonical.message;

/**
 * Branch-length API for linear Gaussian transitions.
 *
 * <p>This interface captures the mathematical kernel of a single Gaussian branch
 * transition without introducing time-index or time-grid concepts. It is therefore
 * suitable for both regularly sampled time series and tree likelihoods, where the
 * primitive quantity is simply a branch length {@code dt}.
 */
public interface GaussianBranchTransitionKernel {

    int getStateDimension();

    void getInitialMean(double[] out);

    void getInitialCovariance(double[][] out);

    void fillTransitionMatrix(double dt, double[][] out);

    void fillTransitionOffset(double dt, double[] out);

    void fillTransitionCovariance(double dt, double[][] out);

    void accumulateSelectionGradient(double dt,
                                     double[][] dLogL_dF,
                                     double[] dLogL_df,
                                     double[] gradientAccumulator);

    default void accumulateSelectionGradientFromCovariance(double dt,
                                                           double[][] dLogL_dV,
                                                           double[] gradientAccumulator) {
        // no-op by default
    }

    default void accumulateDiffusionGradient(double dt,
                                             double[][] dLogL_dV,
                                             double[] gradientAccumulator) {
        // no-op by default
    }
}
