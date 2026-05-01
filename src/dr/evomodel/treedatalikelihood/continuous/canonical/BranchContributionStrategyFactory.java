package dr.evomodel.treedatalikelihood.continuous.canonical;

import dr.evolution.tree.Tree;

final class BranchContributionStrategyFactory {

    private final Tree tree;
    private final int dimension;
    private final CanonicalTreeStateStore stateStore;

    private final BranchContributionStrategy fixedParentObservedTip =
            new FixedParentObservedTipContribution();
    private final BranchContributionStrategy fixedParentPartiallyObservedTip =
            new FixedParentPartiallyObservedTipContribution();
    private final BranchContributionStrategy fixedParentInternalNode =
            new FixedParentInternalNodeContribution();
    private final BranchContributionStrategy observedTip =
            new ObservedTipContribution();
    private final BranchContributionStrategy partiallyObservedTip =
            new PartiallyObservedTipContribution();
    private final BranchContributionStrategy internalNode =
            new InternalNodeContribution();
    private final BranchContributionStrategy missingTrait =
            new MissingTraitContribution();

    BranchContributionStrategyFactory(final Tree tree,
                                      final int dimension,
                                      final CanonicalTreeStateStore stateStore) {
        this.tree = tree;
        this.dimension = dimension;
        this.stateStore = stateStore;
    }

    BranchContributionStrategy select(final int childIndex) {
        final boolean fixedParent =
                stateStore.hasFixedRootValue && tree.isRoot(tree.getParent(tree.getNode(childIndex)));
        if (tree.isExternal(tree.getNode(childIndex))) {
            final CanonicalTipObservation tipObservation = stateStore.tipObservations[childIndex];
            if (tipObservation.observedCount == 0) {
                return missingTrait;
            }
            if (tipObservation.observedCount == dimension) {
                return fixedParent ? fixedParentObservedTip : observedTip;
            }
            return fixedParent ? fixedParentPartiallyObservedTip : partiallyObservedTip;
        }
        return fixedParent ? fixedParentInternalNode : internalNode;
    }
}
