package dr.evomodel.treedatalikelihood.continuous;

/**
 * Minimal tree-side API for branch Gaussian transitions parameterized directly by
 * branch length.
 *
 * <p>This interface is intended as a bridge layer between the tree recursion code and
 * reusable branch kernels originating from the time-series machinery.
 */
public interface GaussianBranchTransitionProvider {

    int getStateDimension();

    void getInitialMean(double[] out);

    void getInitialCovariance(double[][] out);

    void fillBranchTransitionMatrix(double branchLength, double[][] out);

    void fillBranchTransitionOffset(double branchLength, double[] out);

    void fillBranchTransitionCovariance(double branchLength, double[][] out);
}
