package dr.evomodel.continuous.canonical.contribution;

import dr.evomodel.treedatalikelihood.continuousNew.canonical.message.CanonicalGaussianTransition;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.gradient.BranchGradientMessageBundle;
import dr.evomodel.continuous.canonical.workspace.BranchGradientWorkspace;
import dr.evomodel.treedatalikelihood.continuousNew.observationmodel.IdentityCanonicalTipObservationModel;

final class PartiallyObservedTipContribution implements BranchContributionStrategy {

    @Override
    public boolean contributes() {
        return true;
    }

    @Override
    public boolean fillContribution(final BranchGradientMessageBundle messages,
                                    final CanonicalGaussianTransition transition,
                                    final CanonicalBranchContributionAssembler assembler,
                                    final BranchGradientWorkspace workspace) {
        final IdentityCanonicalTipObservationModel tipObservation =
                assembler.identityTipObservationModel(messages);
        assembler.fillContributionForPartiallyObservedTip(
                messages.getAbovePreOrderAtParent(),
                transition,
                tipObservation,
                workspace);
        return true;
    }
}
