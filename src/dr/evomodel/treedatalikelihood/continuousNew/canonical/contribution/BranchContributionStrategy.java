package dr.evomodel.treedatalikelihood.continuousNew.canonical.contribution;

import dr.evomodel.treedatalikelihood.continuousNew.canonical.message.CanonicalGaussianTransition;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.gradient.BranchGradientMessageBundle;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.workspace.BranchGradientWorkspace;

interface BranchContributionStrategy {

    boolean contributes();

    boolean fillContribution(BranchGradientMessageBundle messages,
                             CanonicalGaussianTransition transition,
                             CanonicalBranchContributionAssembler assembler,
                             BranchGradientWorkspace workspace);
}
