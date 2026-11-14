package dr.evomodel.continuous.canonical.contribution;

import dr.evomodel.treedatalikelihood.continuousNew.canonical.gradient.BranchGradientMessageBundle;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.message.CanonicalGaussianTransition;
import dr.evomodel.continuous.canonical.workspace.BranchGradientWorkspace;
import dr.evomodel.treedatalikelihood.continuousNew.observationmodel.CanonicalTipObservationModel;

final class FixedParentGaussianLinkTipBranchContribution implements BranchContributionStrategy {

    @Override
    public boolean contributes() {
        return true;
    }

    @Override
    public boolean fillContribution(final BranchGradientMessageBundle messages,
                                    final CanonicalGaussianTransition transition,
                                    final CanonicalBranchContributionAssembler assembler,
                                    final BranchGradientWorkspace workspace) {
        final CanonicalTipObservationModel observationModel = messages.getTipObservationModel();
        observationModel.fillChildCanonicalState(workspace.state, workspace.observationWorkspace);
        assembler.fillContributionForFixedParentInternalNode(
                transition,
                workspace.state,
                messages.getFixedRootValue(),
                workspace);
        return true;
    }
}
