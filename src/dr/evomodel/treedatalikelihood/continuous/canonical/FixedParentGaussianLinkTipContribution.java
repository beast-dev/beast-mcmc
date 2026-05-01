package dr.evomodel.treedatalikelihood.continuous.canonical;

import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianTransition;
import dr.evomodel.treedatalikelihood.continuous.observationmodel.CanonicalTipObservationModel;

final class FixedParentGaussianLinkTipContribution implements BranchContributionStrategy {

    @Override
    public boolean contributes() {
        return true;
    }

    @Override
    public boolean fillContribution(final int childIndex,
                                    final CanonicalGaussianTransition transition,
                                    final CanonicalBranchContributionAssembler assembler,
                                    final BranchGradientWorkspace workspace) {
        final CanonicalTipObservationModel observationModel =
                assembler.stateStore.tipObservationModels[childIndex];
        observationModel.fillChildCanonicalState(workspace.state, workspace.observationWorkspace);
        assembler.fillContributionForFixedParentInternalNode(transition, workspace.state, workspace);
        return true;
    }
}
