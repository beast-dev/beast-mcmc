package dr.evomodel.continuous.canonical.contribution;

import dr.evomodel.treedatalikelihood.continuousNew.canonical.gradient.BranchGradientMessageBundle;
import dr.evomodel.treedatalikelihood.continuousNew.observationmodel.CanonicalTipObservationModel;
import dr.evomodel.treedatalikelihood.continuousNew.observationmodel.TipObservationMode;

final class BranchContributionStrategyFactory {

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

    BranchContributionStrategyFactory() {
    }

    BranchContributionStrategy select(final BranchGradientMessageBundle messages) {
        final boolean fixedParent = messages.isFixedParent();
        if (messages.isExternalChild()) {
            final CanonicalTipObservationModel observationModel = messages.getTipObservationModel();
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
