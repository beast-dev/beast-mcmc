package dr.evomodel.treedatalikelihood.continuous;

import dr.evomodel.treedatalikelihood.preorder.BranchSufficientStatistics;

interface CanonicalSelectionPullback {

    double[] gradientForBranch(double branchLength,
                               double[] optimum,
                               BranchSufficientStatistics statistics,
                               DenseSelectionGradientProvider denseGradientProvider);

    double[] projectDenseGradient(double[] denseGradient);

    interface DenseSelectionGradientProvider {
        double[] gradient(double branchLength,
                          double[] optimum,
                          BranchSufficientStatistics statistics);
    }
}
