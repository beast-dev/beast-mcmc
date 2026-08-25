package dr.evomodel.branchratemodel;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.TreeParameterModel;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

/**
 * Reward-mixture branch rates driven by one embedded categorical state per
 * branch, exactly like {@link RewardsAwareCategoricalMixtureBranchRates},
 * except atomic states are ordered along the embedded axis by reward value
 * and the continuous category dynamically sits between whichever two atomic
 * states currently bracket that branch's live total.rewards.cts value
 * (see REWARD_CATEGORY_DYNAMIC_ORDERING_PLAN.md). New, opt-in sibling class:
 * does not modify RewardsAwareCategoricalMixtureBranchRates.
 *
 * refreshEmbedding() is deliberately NOT called reactively when
 * total.rewards.cts changes -- an HMC trajectory probes many candidate cts
 * values per operator call, and re-deriving every branch's insertion rank on
 * each of those would reintroduce, per HMC step, the same O(branchCount)
 * redundant-recomputation cost fixed in
 * RewardMixtureCategoricalDiscontinuousPotentialProvider earlier this
 * session. Consumers must call refreshEmbedding() explicitly once per
 * operator call instead, matching where refreshLikelihoodMessages()/
 * refreshEmbedding() are already called today in
 * RewardMixtureCategoricalDiscontinuousPotentialProvider and
 * RewardsMixtureBranchWeightProvider.
 */
public final class RewardsAwareCategoricalMixtureBranchRatesDynamic extends ArbitraryBranchRates
        implements RewardMixtureCategoricalBranchRateModel {

    public static final String ID = "rewardsAwareCategoricalMixtureBranchRatesDynamic";

    private final PerBranchRewardMixtureCategoryDecoder categoryDecoder;
    private final RewardRates rewardRates;
    private final int[] decodedCategoryByParameterIndex;

    public RewardsAwareCategoricalMixtureBranchRatesDynamic(final TreeModel tree,
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
        this.categoryDecoder = new PerBranchRewardMixtureCategoryDecoder(
                categoryParameter,
                categoryCuts,
                ctsParameter,
                rewardRates,
                rewardRates.getStateIndices().getDimension(),
                ctsParameter.getDimension());
        this.decodedCategoryByParameterIndex = new int[ctsParameter.getDimension()];
        refreshCachedDecodedCategories();

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

    public PerBranchRewardMixtureCategoryDecoder getCategoryDecoder() {
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
        return rewardRates.getRawReward(stateIndex);
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
            if (refreshCachedDecodedCategories()) {
                fireModelChanged();
            }
        } else if (variable == categoryDecoder.getCategoryParameter()) {
            if (refreshCachedDecodedCategory(index)) {
                if (index >= 0) {
                    fireModelChanged(variable, getNodeNumberFromParameterIndex(index));
                } else {
                    fireModelChanged();
                }
            }
        } else if (variable == rewardRates.getValues() ||
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
        refreshCachedDecodedCategories();
        super.restoreState();
    }

    private boolean refreshCachedDecodedCategory(final int parameterIndex) {
        if (parameterIndex < 0 || parameterIndex >= decodedCategoryByParameterIndex.length) {
            return refreshCachedDecodedCategories();
        }
        final int currentCategory = categoryDecoder.getCategoryForParameterIndex(parameterIndex);
        final boolean changed = decodedCategoryByParameterIndex[parameterIndex] != currentCategory;
        decodedCategoryByParameterIndex[parameterIndex] = currentCategory;
        return changed;
    }

    private boolean refreshCachedDecodedCategories() {
        boolean changed = false;
        for (int i = 0; i < decodedCategoryByParameterIndex.length; i++) {
            final int currentCategory = categoryDecoder.getCategoryForParameterIndex(i);
            if (decodedCategoryByParameterIndex[i] != currentCategory) {
                changed = true;
                decodedCategoryByParameterIndex[i] = currentCategory;
            }
        }
        return changed;
    }
}
