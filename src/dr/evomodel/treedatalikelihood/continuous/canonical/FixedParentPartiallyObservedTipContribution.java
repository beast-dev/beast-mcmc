package dr.evomodel.treedatalikelihood.continuous.canonical;

import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianTransition;

final class FixedParentPartiallyObservedTipContribution implements BranchContributionStrategy {

    @Override
    public boolean contributes() {
        return true;
    }

    @Override
    public boolean fillContribution(final int childIndex,
                                    final CanonicalGaussianTransition transition,
                                    final CanonicalBranchContributionAssembler assembler,
                                    final BranchGradientWorkspace workspace) {
        final CanonicalTipObservation tipObservation = assembler.stateStore.tipObservations[childIndex];
        final int observedCount = assembler.collectObservationPartition(tipObservation, workspace);
        assembler.fillContributionForFixedParentPartiallyObservedTip(
                transition,
                tipObservation,
                observedCount,
                workspace);
        return true;
    }
}
