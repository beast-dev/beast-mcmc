package dr.evomodel.treedatalikelihood.continuous.canonical.contribution;

import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalBranchMessageContributionUtils;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianMessageOps;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianTransition;
import dr.evomodel.treedatalikelihood.continuous.canonical.workspace.BranchGradientWorkspace;
import dr.evomodel.treedatalikelihood.continuous.observationmodel.CanonicalTipObservationModel;

final class GaussianLinkTipBranchContribution implements BranchContributionStrategy {

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
        CanonicalGaussianMessageOps.buildPairPosterior(
                assembler.stateStore.branchAboveParent[childIndex],
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
