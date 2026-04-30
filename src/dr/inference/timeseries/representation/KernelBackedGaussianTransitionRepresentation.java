package dr.inference.timeseries.representation;

import dr.evomodel.treedatalikelihood.continuous.gaussian.GaussianBranchTransitionKernel;
import dr.inference.timeseries.core.TimeGrid;

/**
 * Thin time-grid adapter over a branch-length Gaussian transition kernel.
 */
public final class KernelBackedGaussianTransitionRepresentation
        implements GaussianTransitionRepresentation {

    private final GaussianBranchTransitionKernel kernel;

    public KernelBackedGaussianTransitionRepresentation(final GaussianBranchTransitionKernel kernel) {
        if (kernel == null) {
            throw new IllegalArgumentException("kernel must not be null");
        }
        this.kernel = kernel;
    }

    @Override
    public int getStateDimension() {
        return kernel.getStateDimension();
    }

    @Override
    public void getInitialMean(final double[] out) {
        kernel.getInitialMean(out);
    }

    @Override
    public void getInitialCovariance(final double[][] out) {
        kernel.getInitialCovariance(out);
    }

    @Override
    public void getTransitionMatrix(final int fromIndex,
                                    final int toIndex,
                                    final TimeGrid timeGrid,
                                    final double[][] out) {
        kernel.fillTransitionMatrix(validatedDelta(timeGrid, fromIndex, toIndex), out);
    }

    @Override
    public void getTransitionOffset(final int fromIndex,
                                    final int toIndex,
                                    final TimeGrid timeGrid,
                                    final double[] out) {
        kernel.fillTransitionOffset(validatedDelta(timeGrid, fromIndex, toIndex), out);
    }

    @Override
    public void getTransitionCovariance(final int fromIndex,
                                        final int toIndex,
                                        final TimeGrid timeGrid,
                                        final double[][] out) {
        kernel.fillTransitionCovariance(validatedDelta(timeGrid, fromIndex, toIndex), out);
    }

    @Override
    public void accumulateSelectionGradient(final int fromIndex,
                                            final int toIndex,
                                            final TimeGrid timeGrid,
                                            final double[][] dLogL_dF,
                                            final double[] dLogL_df,
                                            final double[] gradientAccumulator) {
        kernel.accumulateSelectionGradient(
                validatedDelta(timeGrid, fromIndex, toIndex),
                dLogL_dF,
                dLogL_df,
                gradientAccumulator);
    }

    @Override
    public void accumulateSelectionGradientFromCovariance(final int fromIndex,
                                                          final int toIndex,
                                                          final TimeGrid timeGrid,
                                                          final double[][] dLogL_dV,
                                                          final double[] gradientAccumulator) {
        kernel.accumulateSelectionGradientFromCovariance(
                validatedDelta(timeGrid, fromIndex, toIndex),
                dLogL_dV,
                gradientAccumulator);
    }

    @Override
    public void accumulateDiffusionGradient(final int fromIndex,
                                            final int toIndex,
                                            final TimeGrid timeGrid,
                                            final double[][] dLogL_dV,
                                            final double[] gradientAccumulator) {
        kernel.accumulateDiffusionGradient(
                validatedDelta(timeGrid, fromIndex, toIndex),
                dLogL_dV,
                gradientAccumulator);
    }

    private static double validatedDelta(final TimeGrid timeGrid, final int fromIndex, final int toIndex) {
        final double dt = timeGrid.getDelta(fromIndex, toIndex);
        if (!(dt > 0.0)) {
            throw new IllegalArgumentException("Time increments must be strictly positive");
        }
        return dt;
    }
}
