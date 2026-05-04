package dr.evomodel.treedatalikelihood.continuous.canonical.contribution;

import dr.evolution.tree.Tree;
import dr.evomodel.treedatalikelihood.continuous.canonical.traversal.CanonicalTreeStateStore;
import dr.evomodel.treedatalikelihood.continuous.observationmodel.CanonicalTipObservationModel;
import dr.evomodel.treedatalikelihood.continuous.observationmodel.TipObservationMode;

final class BranchContributionStrategyFactory {

    private final Tree tree;
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
    private final BranchContributionStrategy gaussianLinkTip =
            new GaussianLinkTipBranchContribution();
    private final BranchContributionStrategy fixedParentGaussianLinkTip =
            new FixedParentGaussianLinkTipBranchContribution();

    BranchContributionStrategyFactory(final Tree tree,
                                      final int dimension,
                                      final CanonicalTreeStateStore stateStore) {
        this.tree = tree;
        this.stateStore = stateStore;
    }

    BranchContributionStrategy select(final int childIndex) {
        final boolean fixedParent =
                stateStore.hasFixedRootValue && tree.isRoot(tree.getParent(tree.getNode(childIndex)));
        if (tree.isExternal(tree.getNode(childIndex))) {
            final CanonicalTipObservationModel observationModel =
                    stateStore.tipObservationModels[childIndex];
            if (observationModel.getMode() == TipObservationMode.GAUSSIAN_LINK) {
                return fixedParent ? fixedParentGaussianLinkTip : gaussianLinkTip;
            }
            switch (observationModel.getMode()) {
                case MISSING:
                    return missingTrait;
                case EXACT_IDENTITY:
                    return fixedParent ? fixedParentObservedTip : observedTip;
                case PARTIAL_EXACT_IDENTITY:
                    return fixedParent ? fixedParentPartiallyObservedTip : partiallyObservedTip;
                default:
                    throw new IllegalStateException("Unsupported tip observation mode: "
                            + observationModel.getMode());
            }
        }
        return fixedParent ? fixedParentInternalNode : internalNode;
    }
}
