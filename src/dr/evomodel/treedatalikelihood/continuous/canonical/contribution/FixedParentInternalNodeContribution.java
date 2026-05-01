package dr.evomodel.treedatalikelihood.continuous.canonical.contribution;

import dr.evomodel.treedatalikelihood.continuous.canonical.*;
import dr.evomodel.treedatalikelihood.continuous.canonical.traversal.*;
import dr.evomodel.treedatalikelihood.continuous.canonical.workspace.*;

import dr.evomodel.treedatalikelihood.continuous.canonical.workspace.*;

import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianTransition;

final class FixedParentInternalNodeContribution implements BranchContributionStrategy {

    @Override
    public boolean contributes() {
        return true;
    }

    @Override
    public boolean fillContribution(final int childIndex,
                                    final CanonicalGaussianTransition transition,
                                    final CanonicalBranchContributionAssembler assembler,
                                    final BranchGradientWorkspace workspace) {
        assembler.fillContributionForFixedParentInternalNode(
                transition,
                assembler.stateStore.postOrder[childIndex],
                workspace);
        return true;
    }
}
