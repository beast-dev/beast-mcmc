package dr.evomodel.branchratemodel;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.TreeParameterModel;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

/**
 * Reward-mixture branch rates driven by one embedded categorical state per
 * branch.
 *
 * Category 0 uses the latent continuous reward parameter. Categories 1..K use
 * reward atom states 0..K-1.
 */
public final class RewardsAwareCategoricalMixtureBranchRates extends ArbitraryBranchRates
        implements RewardMixtureBranchRateModel {

    public static final String ID = "rewardsAwareCategoricalMixtureBranchRates";

    private final RewardMixtureCategoryDecoder categoryDecoder;
    private final RewardRates rewardRates;

    public RewardsAwareCategoricalMixtureBranchRates(final TreeModel tree,
                                                     final Parameter ctsParameter,
                                                     final Parameter categoryParameter,
                                                     final Parameter categoryCuts,
                                                     final RewardRates rewardRates,
                                                     final BranchRateTransform transform,
                                                     final boolean setRates,
                                                     final TreeParameterModel.Type includeRoot) {
        super(ID, tree, ctsParameter,
                transform == null ? new BranchRateTransform.None() : transform,
                setRates, includeRoot);

        if (rewardRates == null) {
            throw new IllegalArgumentException("rewardRates must be non-null");
        }

        this.rewardRates = rewardRates;
        this.categoryDecoder = new RewardMixtureCategoryDecoder(
                categoryParameter,
                categoryCuts,
                rewardRates.getStateIndices().getDimension(),
                ctsParameter.getDimension());

        addVariable(categoryParameter);
        addVariable(categoryCuts);
        addVariable(rewardRates.getValues());
        addVariable(rewardRates.getVaryingValues());
        addVariable(rewardRates.getStateIndices());
    }

    @Override
    public double getUntransformedBranchRate(final Tree tree, final NodeRef node) {
        final int p = getParameterIndexFromNode(node);
        if (categoryDecoder.isAtomic(p)) {
            return getRawRewardForAtomState(categoryDecoder.getAtomicState(p));
        }
        return getContinuousRawReward(tree, node);
    }

    @Override
    public double getBranchRateDifferential(final Tree tree, final NodeRef node) {
        final int p = getParameterIndexFromNode(node);
        return categoryDecoder.isAtomic(p) ? 0.0 : super.getBranchRateDifferential(tree, node);
    }

    @Override
    public double getBranchRateSecondDifferential(final Tree tree, final NodeRef node) {
        final int p = getParameterIndexFromNode(node);
        return categoryDecoder.isAtomic(p) ? 0.0 : super.getBranchRateSecondDifferential(tree, node);
    }

    @Override
    public double getContinuousRawReward(final Tree tree, final NodeRef node) {
        return super.getUntransformedBranchRate(tree, node);
    }

    @Override
    public double getBranchRateForRawReward(final Tree tree, final NodeRef node, final double rawReward) {
        return getTransform().transform(rawReward, tree, node);
    }

    public RewardMixtureCategoryDecoder getCategoryDecoder() {
        return categoryDecoder;
    }

    public Parameter getCategoryParameter() {
        return categoryDecoder.getCategoryParameter();
    }

    public Parameter getCategoryCutParameter() {
        return categoryDecoder.getCutParameter();
    }

    public RewardRates getRewardRates() {
        return rewardRates;
    }

    public double getRawRewardForAtomState(final int stateIndex) {
        if (stateIndex < 0 || stateIndex >= rewardRates.getStateIndices().getDimension()) {
            throw new IllegalArgumentException("stateIndex out of range: " + stateIndex);
        }
        final int rewardRateIndex = (int) Math.round(rewardRates.getStateIndices().getParameterValue(stateIndex));
        if (rewardRateIndex < 0 || rewardRateIndex >= rewardRates.getValues().getDimension()) {
            throw new IllegalArgumentException(
                    "Reward-rate mapping for state " + stateIndex + " points outside rewardRates: " +
                            rewardRateIndex);
        }
        return rewardRates.getValues().getParameterValue(rewardRateIndex);
    }

    public Parameter getRewardRatesValues() {
        return rewardRates.getValues();
    }

    public Parameter getRewardRatesInternal() {
        return rewardRates.getVaryingValues();
    }

    public Parameter getRewardRatesMapping() {
        return rewardRates.getStateIndices();
    }

    @Override
    protected void handleVariableChangedEvent(final Variable variable,
                                              final int index,
                                              final Parameter.ChangeType type) {
        if (variable == categoryDecoder.getCutParameter()) {
            categoryDecoder.refreshEmbedding();
            fireModelChanged();
        } else if (variable == categoryDecoder.getCategoryParameter() ||
                variable == rewardRates.getValues() ||
                variable == rewardRates.getVaryingValues() ||
                variable == rewardRates.getStateIndices()) {
            fireModelChanged();
        } else {
            super.handleVariableChangedEvent(variable, index, type);
        }
    }

    @Override
    protected void restoreState() {
        categoryDecoder.refreshEmbedding();
        super.restoreState();
    }
}
