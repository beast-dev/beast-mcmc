package dr.evomodel.treedatalikelihood.continuousNew.canonical.contribution;

import dr.evomodel.treedatalikelihood.continuousNew.canonical.message.CanonicalBranchMessageContributionUtils;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.message.CanonicalGaussianMessageOps;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.message.CanonicalGaussianTransition;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.gradient.BranchGradientMessageBundle;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.workspace.BranchGradientWorkspace;

final class InternalNodeContribution implements BranchContributionStrategy {

    @Override
    public boolean contributes() {
        return true;
    }

    @Override
    public boolean fillContribution(final BranchGradientMessageBundle messages,
                                    final CanonicalGaussianTransition transition,
                                    final CanonicalBranchContributionAssembler assembler,
                                    final BranchGradientWorkspace workspace) {
        CanonicalGaussianMessageOps.buildPairPosterior(
                messages.getAbovePreOrderAtParent(),
                transition,
                messages.getBelowPostOrderAtChild(),
                workspace.pairState);
        CanonicalBranchMessageContributionUtils.fillFromPairState(
                workspace.pairState,
                workspace.contributionWorkspace,
                workspace.contribution);
        return true;
    }
}
