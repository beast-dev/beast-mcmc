package dr.evomodel.continuous.canonical.contribution;

import dr.evomodel.treedatalikelihood.continuousNew.canonical.message.CanonicalGaussianTransition;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.gradient.BranchGradientMessageBundle;
import dr.evomodel.continuous.canonical.workspace.BranchGradientWorkspace;

final class FixedParentInternalNodeContribution implements BranchContributionStrategy {

    @Override
    public boolean contributes() {
        return true;
    }

    @Override
    public boolean fillContribution(final BranchGradientMessageBundle messages,
                                    final CanonicalGaussianTransition transition,
                                    final CanonicalBranchContributionAssembler assembler,
                                    final BranchGradientWorkspace workspace) {
        assembler.fillContributionForFixedParentInternalNode(
                transition,
                messages.getBelowPostOrderAtChild(),
                messages.getFixedRootValue(),
                workspace);
        return true;
    }
}
