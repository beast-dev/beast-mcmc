package dr.evomodel.treedatalikelihood.continuous.canonical;

import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalBranchMessageContributionUtils;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianMessageOps;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianTransition;

final class InternalNodeContribution implements BranchContributionStrategy {

    @Override
    public boolean contributes() {
        return true;
    }

    @Override
    public boolean fillContribution(final int childIndex,
                                    final CanonicalGaussianTransition transition,
                                    final CanonicalBranchContributionAssembler assembler,
                                    final BranchGradientWorkspace workspace) {
        CanonicalGaussianMessageOps.buildPairPosterior(
                assembler.stateStore.branchAboveParent[childIndex],
                transition,
                assembler.stateStore.postOrder[childIndex],
                workspace.pairState);
        CanonicalBranchMessageContributionUtils.fillFromPairState(
                workspace.pairState,
                workspace.contributionWorkspace,
                workspace.contribution);
        return true;
    }
}
