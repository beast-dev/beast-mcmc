package dr.evomodel.treedatalikelihood.continuous.integration;

import dr.evomodel.treedatalikelihood.continuous.gaussian.CanonicalGaussianState;
import dr.evomodel.treedatalikelihood.continuous.gaussian.message.CanonicalBranchMessageContribution;
import dr.evomodel.treedatalikelihood.continuous.gaussian.message.CanonicalBranchMessageContributionUtils;
import dr.evomodel.treedatalikelihood.continuous.gaussian.message.CanonicalLocalTransitionAdjoints;
import dr.evomodel.treedatalikelihood.continuous.gaussian.message.CanonicalTransitionAdjointUtils;

final class BranchAdjointWorkspace {
    final CanonicalGaussianState combinedState;
    final CanonicalGaussianState parentPosterior;
    final CanonicalGaussianState pairState;
    final CanonicalBranchMessageContributionUtils.Workspace contributionWorkspace;
    final CanonicalTransitionAdjointUtils.Workspace transitionAdjointWorkspace;
    final CanonicalBranchMessageContribution contribution;
    final CanonicalLocalTransitionAdjoints adjoints;
    final int[] observedIndexScratch;
    final int[] missingIndexScratch;
    final int[] reducedIndexByTraitScratch;
    final double[] mean;
    final double[] mean2;
    final double[][] covariance;
    final double[] varianceFlat;
    final double[] precisionFlat;
    final double[] reducedPrecisionFlatScratch;
    final double[] reducedCovarianceFlatScratch;
    final double[] reducedInformationScratch;
    final double[] reducedMeanScratch;

    BranchAdjointWorkspace(final int dim) {
        this.combinedState = new CanonicalGaussianState(dim);
        this.parentPosterior = new CanonicalGaussianState(dim);
        this.pairState = new CanonicalGaussianState(2 * dim);
        this.contributionWorkspace = new CanonicalBranchMessageContributionUtils.Workspace(dim);
        this.transitionAdjointWorkspace = new CanonicalTransitionAdjointUtils.Workspace(dim);
        this.contribution = new CanonicalBranchMessageContribution(dim);
        this.adjoints = new CanonicalLocalTransitionAdjoints(dim);
        this.observedIndexScratch = new int[dim];
        this.missingIndexScratch = new int[dim];
        this.reducedIndexByTraitScratch = new int[dim];
        this.mean = new double[dim];
        this.mean2 = new double[dim];
        this.covariance = new double[dim][dim];
        this.varianceFlat = new double[dim * dim];
        this.precisionFlat = new double[dim * dim];
        this.reducedPrecisionFlatScratch = new double[4 * dim * dim];
        this.reducedCovarianceFlatScratch = new double[4 * dim * dim];
        this.reducedInformationScratch = new double[2 * dim];
        this.reducedMeanScratch = new double[2 * dim];
    }
}
