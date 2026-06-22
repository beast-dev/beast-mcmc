package dr.evomodel.treedatalikelihood.continuous.canonical.contribution;

import dr.evomodel.treedatalikelihood.continuous.canonical.gradient.BranchGradientMessageBundle;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianTransition;
import dr.evomodel.treedatalikelihood.continuous.canonical.workspace.BranchGradientWorkspace;
import dr.evomodel.treedatalikelihood.continuous.observationmodel.CanonicalTipObservationModel;

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
