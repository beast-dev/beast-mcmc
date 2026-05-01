package dr.evomodel.treedatalikelihood.continuous.canonical.contribution;

import dr.evomodel.treedatalikelihood.continuous.canonical.*;
import dr.evomodel.treedatalikelihood.continuous.canonical.traversal.*;
import dr.evomodel.treedatalikelihood.continuous.canonical.workspace.*;

import dr.evomodel.treedatalikelihood.continuous.canonical.workspace.*;

import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianTransition;
import dr.evomodel.treedatalikelihood.continuous.observationmodel.IdentityCanonicalTipObservationModel;

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
        final IdentityCanonicalTipObservationModel tipObservation = assembler.identityTipObservationModel(childIndex);
        assembler.fillContributionForFixedParentPartiallyObservedTip(
                transition,
                tipObservation,
                workspace);
        return true;
    }
}
