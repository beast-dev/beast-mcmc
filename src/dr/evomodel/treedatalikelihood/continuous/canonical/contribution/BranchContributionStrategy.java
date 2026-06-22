package dr.evomodel.treedatalikelihood.continuous.canonical.contribution;

import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianTransition;
import dr.evomodel.treedatalikelihood.continuous.canonical.gradient.BranchGradientMessageBundle;
import dr.evomodel.treedatalikelihood.continuous.canonical.workspace.BranchGradientWorkspace;

interface BranchContributionStrategy {

    boolean contributes();

    boolean fillContribution(BranchGradientMessageBundle messages,
                             CanonicalGaussianTransition transition,
                             CanonicalBranchContributionAssembler assembler,
                             BranchGradientWorkspace workspace);
}
