package dr.evomodel.treedatalikelihood.continuousNew.canonical.contribution;

import dr.evomodel.treedatalikelihood.continuousNew.canonical.gradient.BranchGradientMessageBundle;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.message.CanonicalBranchMessageContributionUtils;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.message.CanonicalGaussianMessageOps;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.message.CanonicalGaussianTransition;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.workspace.BranchGradientWorkspace;
import dr.evomodel.treedatalikelihood.continuousNew.observationmodel.CanonicalTipObservationModel;

final class GaussianLinkTipBranchContribution implements BranchContributionStrategy {

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
        CanonicalGaussianMessageOps.buildPairPosterior(
                messages.getAbovePreOrderAtParent(),
                transition,
                workspace.state,
                workspace.pairState);
        CanonicalBranchMessageContributionUtils.fillFromPairState(
                workspace.pairState,
                workspace.contributionWorkspace,
                workspace.contribution);
        return true;
    }
}
