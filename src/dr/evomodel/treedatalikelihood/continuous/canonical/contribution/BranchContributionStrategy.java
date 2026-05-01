package dr.evomodel.treedatalikelihood.continuous.canonical.contribution;

import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianTransition;
import dr.evomodel.treedatalikelihood.continuous.canonical.workspace.BranchGradientWorkspace;

interface BranchContributionStrategy {

    boolean contributes();

    boolean fillContribution(int childIndex,
                             CanonicalGaussianTransition transition,
                             CanonicalBranchContributionAssembler assembler,
                             BranchGradientWorkspace workspace);
}
