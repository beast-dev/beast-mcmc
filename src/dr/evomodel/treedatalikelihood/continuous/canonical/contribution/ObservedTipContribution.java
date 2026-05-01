package dr.evomodel.treedatalikelihood.continuous.canonical.contribution;

import dr.evomodel.treedatalikelihood.continuous.canonical.*;
import dr.evomodel.treedatalikelihood.continuous.canonical.traversal.*;
import dr.evomodel.treedatalikelihood.continuous.canonical.workspace.*;

import dr.evomodel.treedatalikelihood.continuous.canonical.workspace.*;

import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianTransition;

final class ObservedTipContribution implements BranchContributionStrategy {

    @Override
    public boolean contributes() {
        return true;
    }

    @Override
    public boolean fillContribution(final int childIndex,
                                    final CanonicalGaussianTransition transition,
                                    final CanonicalBranchContributionAssembler assembler,
                                    final BranchGradientWorkspace workspace) {
        assembler.fillContributionForObservedTip(
                assembler.stateStore.branchAboveParent[childIndex],
                transition,
                assembler.identityTipObservationModel(childIndex),
                workspace);
        return true;
    }
}
