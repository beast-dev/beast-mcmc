package dr.evomodel.treedatalikelihood.continuousNew.canonical.contribution;

import dr.evomodel.treedatalikelihood.continuousNew.canonical.gradient.BranchGradientMessageBundle;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.message.CanonicalGaussianTransition;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.workspace.BranchGradientWorkspace;

final class MissingTraitContribution implements BranchContributionStrategy {

    @Override
    public boolean contributes() {
        return false;
    }

    @Override
    public boolean fillContribution(final BranchGradientMessageBundle messages,
                                    final CanonicalGaussianTransition transition,
                                    final CanonicalBranchContributionAssembler assembler,
                                    final BranchGradientWorkspace workspace) {
        return false;
    }
}
