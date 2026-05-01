package dr.evomodel.treedatalikelihood.continuous.canonical;

import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianTransition;

interface BranchContributionStrategy {

    boolean contributes();

    boolean fillContribution(int childIndex,
                             CanonicalGaussianTransition transition,
                             CanonicalBranchContributionAssembler assembler,
                             BranchGradientWorkspace workspace);
}
