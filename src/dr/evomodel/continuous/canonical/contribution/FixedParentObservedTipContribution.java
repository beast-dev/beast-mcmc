package dr.evomodel.continuous.canonical.contribution;

import dr.evomodel.treedatalikelihood.continuousNew.canonical.message.CanonicalGaussianTransition;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.gradient.BranchGradientMessageBundle;
import dr.evomodel.continuous.canonical.workspace.BranchGradientWorkspace;

final class FixedParentObservedTipContribution implements BranchContributionStrategy {

    @Override
    public boolean contributes() {
        return true;
    }

    @Override
    public boolean fillContribution(final BranchGradientMessageBundle messages,
                                    final CanonicalGaussianTransition transition,
                                    final CanonicalBranchContributionAssembler assembler,
                                    final BranchGradientWorkspace workspace) {
        assembler.fillContributionForFixedParentObservedTip(
                transition,
                assembler.identityTipObservationModel(messages),
                messages.getFixedRootValue(),
                workspace);
        return true;
    }
}
