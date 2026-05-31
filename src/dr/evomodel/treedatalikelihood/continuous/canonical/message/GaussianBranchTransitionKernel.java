package dr.evomodel.treedatalikelihood.continuous.canonical.message;

import dr.evomodel.treedatalikelihood.continuous.canonical.math.MatrixOps;

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

    default void getInitialCovarianceFlat(final double[] out) {
        final int dim = getStateDimension();
        final double[][] matrix = new double[dim][dim];
        getInitialCovariance(matrix);
        MatrixOps.toFlat(matrix, out, dim);
    }

    void fillTransitionMatrix(double dt, double[][] out);

    default void fillTransitionMatrixFlat(final double dt,
                                          final double[] out) {
        final int dim = getStateDimension();
        final double[][] matrix = new double[dim][dim];
        fillTransitionMatrix(dt, matrix);
        MatrixOps.toFlat(matrix, out, dim);
    }

    void fillTransitionOffset(double dt, double[] out);

    void fillTransitionCovariance(double dt, double[][] out);

    default void fillTransitionCovarianceFlat(final double dt,
                                              final double[] out) {
        final int dim = getStateDimension();
        final double[][] matrix = new double[dim][dim];
        fillTransitionCovariance(dt, matrix);
        MatrixOps.toFlat(matrix, out, dim);
    }

    void accumulateSelectionGradient(double dt,
                                     double[][] dLogL_dF,
                                     double[] dLogL_df,
                                     double[] gradientAccumulator);

    default void accumulateSelectionGradientFlat(final double dt,
                                                 final double[] dLogL_dF,
                                                 final double[] dLogL_df,
                                                 final double[] gradientAccumulator) {
        final int dim = getStateDimension();
        final double[][] matrix = new double[dim][dim];
        MatrixOps.fromFlat(dLogL_dF, matrix, dim);
        accumulateSelectionGradient(dt, matrix, dLogL_df, gradientAccumulator);
    }

    default void accumulateSelectionGradientFromCovariance(double dt,
                                                           double[][] dLogL_dV,
                                                           double[] gradientAccumulator) {
        // no-op by default
    }

    default void accumulateSelectionGradientFromCovarianceFlat(final double dt,
                                                              final double[] dLogL_dV,
                                                              final double[] gradientAccumulator) {
        final int dim = getStateDimension();
        final double[][] matrix = new double[dim][dim];
        MatrixOps.fromFlat(dLogL_dV, matrix, dim);
        accumulateSelectionGradientFromCovariance(dt, matrix, gradientAccumulator);
    }

    default void accumulateDiffusionGradient(double dt,
                                             double[][] dLogL_dV,
                                             double[] gradientAccumulator) {
        // no-op by default
    }

    default void accumulateDiffusionGradientFlat(final double dt,
                                                 final double[] dLogL_dV,
                                                 final double[] gradientAccumulator) {
        final int dim = getStateDimension();
        final double[][] matrix = new double[dim][dim];
        MatrixOps.fromFlat(dLogL_dV, matrix, dim);
        accumulateDiffusionGradient(dt, matrix, gradientAccumulator);
    }
}
