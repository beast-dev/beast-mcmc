package dr.inference.timeseries.representation;

import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianTransition;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.GaussianBranchTransitionKernel;
import dr.inference.timeseries.core.TimeGrid;

/**
 * Thin time-grid adapter over a branch-length Gaussian transition kernel.
 */
public final class KernelBackedGaussianTransitionRepresentation
        implements CachedGaussianTransitionRepresentation {

    private final GaussianBranchTransitionKernel kernel;
    private final RepeatedDeltaCanonicalTransitionCache repeatedDeltaCache;

    public KernelBackedGaussianTransitionRepresentation(final GaussianBranchTransitionKernel kernel) {
        this(kernel, null);
    }

    public KernelBackedGaussianTransitionRepresentation(final OUProcessModel processModel) {
        this(processModel, processModel);
    }

    private KernelBackedGaussianTransitionRepresentation(final GaussianBranchTransitionKernel kernel,
                                                         final OUProcessModel processModel) {
        if (kernel == null) {
            throw new IllegalArgumentException("kernel must not be null");
        }
        this.kernel = kernel;
        this.repeatedDeltaCache = new RepeatedDeltaCanonicalTransitionCache(kernel, processModel);
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
        repeatedDeltaCache.fillTransitionMatrix(validatedDelta(timeGrid, fromIndex, toIndex), out);
    }

    @Override
    public void getTransitionOffset(final int fromIndex,
                                    final int toIndex,
                                    final TimeGrid timeGrid,
                                    final double[] out) {
        repeatedDeltaCache.fillTransitionOffset(validatedDelta(timeGrid, fromIndex, toIndex), out);
    }

    @Override
    public void getTransitionCovariance(final int fromIndex,
                                        final int toIndex,
                                        final TimeGrid timeGrid,
                                        final double[][] out) {
        repeatedDeltaCache.fillTransitionCovariance(validatedDelta(timeGrid, fromIndex, toIndex), out);
    }

    @Override
    public void prepareTimeGrid(final TimeGrid timeGrid) {
        repeatedDeltaCache.prepareTimeGrid(timeGrid);
    }

    @Override
    public void getCanonicalTransition(final int fromIndex,
                                       final int toIndex,
                                       final TimeGrid timeGrid,
                                       final CanonicalGaussianTransition out) {
        repeatedDeltaCache.fillCanonicalTransition(validatedDelta(timeGrid, fromIndex, toIndex), out);
    }

    @Override
    public RepeatedDeltaCacheStatistics getCacheStatistics() {
        return repeatedDeltaCache.getStatistics();
    }

    @Override
    public void makeDirty() {
        repeatedDeltaCache.makeDirty();
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
        return RepeatedDeltaCanonicalTransitionCache.validatedDelta(timeGrid, fromIndex, toIndex);
    }
}
