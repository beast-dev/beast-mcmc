package dr.evomodel.treedatalikelihood.continuous.canonical.contribution;

import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianTransition;
import dr.evomodel.treedatalikelihood.continuous.canonical.workspace.BranchGradientWorkspace;
import dr.evomodel.treedatalikelihood.continuous.observationmodel.IdentityCanonicalTipObservationModel;

final class PartiallyObservedTipContribution implements BranchContributionStrategy {

    @Override
    public boolean contributes() {
        return true;
    }

    @Override
    public boolean fillContribution(final int childIndex,
                                    final CanonicalGaussianTransition transition,
                                    final CanonicalBranchContributionAssembler assembler,
                                    final BranchGradientWorkspace workspace) {
        final IdentityCanonicalTipObservationModel tipObservation = assembler.identityTipObservationModel(childIndex);
        assembler.fillContributionForPartiallyObservedTip(
                assembler.stateStore.branchAboveParent[childIndex],
                transition,
                tipObservation,
                workspace);
        return true;
    }
}
