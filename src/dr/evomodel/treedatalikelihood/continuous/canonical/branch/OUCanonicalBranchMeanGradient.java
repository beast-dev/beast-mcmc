package dr.evomodel.treedatalikelihood.continuous.canonical.branch;

import dr.evomodel.treedatalikelihood.continuous.CanonicalDebugOptions;
import dr.evomodel.treedatalikelihood.continuous.CanonicalGradientFallbackPolicy;
import dr.evomodel.treedatalikelihood.continuous.OUGaussianBranchTransitionProvider;

import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalLocalTransitionAdjoints;
import dr.evomodel.treedatalikelihood.preorder.BranchSufficientStatistics;

final class OUCanonicalBranchMeanGradient {

    private final OUCanonicalBranchWiring branchWiring;
    private final OUCanonicalBranchFiniteDifference finiteDifference;
    private final OUCanonicalBranchGradientScratch scratch;

    OUCanonicalBranchMeanGradient(final OUCanonicalBranchWiring branchWiring,
                                  final OUCanonicalBranchFiniteDifference finiteDifference,
                                  final OUCanonicalBranchGradientScratch scratch) {
        this.branchWiring = branchWiring;
        this.finiteDifference = finiteDifference;
        this.scratch = scratch;
    }

    double[] getGradientWrtBranchDrift(final BranchSufficientStatistics statistics) {
        final CanonicalLocalTransitionAdjoints localAdjoints = scratch.localAdjoints;
        branchWiring.fillLocalAdjoints(statistics, localAdjoints);
        branchWiring.accumulateBranchMeanGradient(
                statistics.getBranch(),
                localAdjoints.dLogL_df,
                scratch.branchMeanGradient);
        return scratch.branchMeanGradient.clone();
    }

    double[] getGradientWrtBranchDrift(final double branchLength,
                                       final double[] optimum,
                                       final BranchSufficientStatistics statistics,
                                       final int nodeNumber) {
        final CanonicalLocalTransitionAdjoints localAdjoints = scratch.localAdjoints;
        branchWiring.fillLocalAdjoints(statistics, localAdjoints);
        branchWiring.accumulateBranchMeanGradient(
                statistics.getBranch(),
                localAdjoints.dLogL_df,
                scratch.branchMeanGradient);
        finiteDifference.maybeEmitLocalMeanFDDebug(
                branchLength,
                optimum,
                statistics,
                nodeNumber,
                scratch.branchMeanGradient,
                localAdjoints);
        return scratch.branchMeanGradient.clone();
    }
}
