package dr.evomodel.treedatalikelihood.continuous.canonical.message;

import dr.evomodel.treedatalikelihood.continuous.canonical.math.MatrixOps;

/**
 * Branch-length API for linear Gaussian transitions.
 *
 * <p>This interface captures the mathematical kernel of a single Gaussian branch
 * transition without introducing time-index or time-grid concepts. It is therefore
 * suitable for both regularly sampled time series and tree likelihoods, where the
 * primitive quantity is simply a branch length {@code dt}.
 *
 * <p>The primary methods use row-major flattened matrices. The {@code double[][]}
 * methods are boundary adapters for older tree-side callers.
 */
public interface GaussianBranchTransitionKernel {

    int getStateDimension();

    void getInitialMean(double[] out);

    void getInitialCovarianceFlat(double[] out);

    default void getInitialCovariance(final double[][] out) {
        final int dim = getStateDimension();
        final double[] flat = new double[dim * dim];
        getInitialCovarianceFlat(flat);
        MatrixOps.fromFlat(flat, out, dim);
    }

    void fillTransitionMatrixFlat(double dt, double[] out);

    default void fillTransitionMatrix(final double dt, final double[][] out) {
        final int dim = getStateDimension();
        final double[] flat = new double[dim * dim];
        fillTransitionMatrixFlat(dt, flat);
        MatrixOps.fromFlat(flat, out, dim);
    }

    void fillTransitionOffset(double dt, double[] out);

    void fillTransitionCovarianceFlat(double dt, double[] out);

    default void fillTransitionCovariance(final double dt, final double[][] out) {
        final int dim = getStateDimension();
        final double[] flat = new double[dim * dim];
        fillTransitionCovarianceFlat(dt, flat);
        MatrixOps.fromFlat(flat, out, dim);
    }

    void accumulateSelectionGradientFlat(double dt,
                                         double[] dLogL_dF,
                                         double[] dLogL_df,
                                         double[] gradientAccumulator);

    default void accumulateSelectionGradientFromCovarianceFlat(final double dt,
                                                              final double[] dLogL_dV,
                                                              final double[] gradientAccumulator) {
        // no-op by default
    }

    default void accumulateDiffusionGradientFlat(final double dt,
                                                 final double[] dLogL_dV,
                                                 final double[] gradientAccumulator) {
        // no-op by default
    }
}
