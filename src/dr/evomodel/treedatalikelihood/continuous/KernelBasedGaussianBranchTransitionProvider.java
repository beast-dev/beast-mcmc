package dr.evomodel.treedatalikelihood.continuous;

import dr.evomodel.treedatalikelihood.continuous.gaussian.GaussianBranchTransitionKernel;

/**
 * Tree-side adapter over a shared {@link GaussianBranchTransitionKernel}.
 *
 * <p>This class is intentionally thin: it lets the tree continuous-data likelihood
 * consume the same branch-length Gaussian transition machinery used by the
 * time-series code, without introducing fake time grids or indices.
 */
public final class KernelBasedGaussianBranchTransitionProvider
        implements GaussianBranchTransitionProvider {

    private final GaussianBranchTransitionKernel kernel;

    public KernelBasedGaussianBranchTransitionProvider(final GaussianBranchTransitionKernel kernel) {
        if (kernel == null) {
            throw new IllegalArgumentException("kernel must not be null");
        }
        this.kernel = kernel;
    }

    public GaussianBranchTransitionKernel getKernel() {
        return kernel;
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
    public void fillBranchTransitionMatrix(final double branchLength, final double[][] out) {
        validateBranchLength(branchLength);
        kernel.fillTransitionMatrix(branchLength, out);
    }

    @Override
    public void fillBranchTransitionOffset(final double branchLength, final double[] out) {
        validateBranchLength(branchLength);
        kernel.fillTransitionOffset(branchLength, out);
    }

    @Override
    public void fillBranchTransitionCovariance(final double branchLength, final double[][] out) {
        validateBranchLength(branchLength);
        kernel.fillTransitionCovariance(branchLength, out);
    }

    private static void validateBranchLength(final double branchLength) {
        if (!(branchLength > 0.0)) {
            throw new IllegalArgumentException("Branch lengths must be strictly positive");
        }
    }
}
