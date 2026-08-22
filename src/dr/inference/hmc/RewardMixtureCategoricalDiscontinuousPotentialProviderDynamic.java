package dr.inference.hmc;

import dr.evomodel.branchmodel.RewardsAwareBranchModel;
import dr.evomodel.branchratemodel.PerBranchRewardMixtureCategoryDecoder;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.inference.model.Parameter;
import dr.inference.operators.RewardsMixtureBranchResamplingHelper;
import dr.inference.operators.RewardsMixtureBranchWeightProvider;

/**
 * Discontinuous-HMC potential provider for reward-mixture branch categories,
 * exactly like {@link RewardMixtureCategoricalDiscontinuousPotentialProvider}
 * except decoding goes through a branch-specific
 * {@link PerBranchRewardMixtureCategoryDecoder} (atomic states ordered by
 * reward value, continuous category tracks each branch's live
 * total.rewards.cts value) instead of one embedding shared by every branch.
 * New, opt-in sibling class: does not modify the existing provider.
 *
 * Takes a pre-built decoder rather than raw categoryParameter/categoryCuts so
 * it can share the same decoder instance -- and the same
 * once-per-operator-call refreshEmbedding() cost -- as the
 * RewardsAwareCategoricalMixtureBranchRatesDynamic and RewardsAwareBranchModel
 * that already reference it, rather than each maintaining its own.
 */
public final class RewardMixtureCategoricalDiscontinuousPotentialProviderDynamic
        implements DiscontinuousPotentialProvider {

    private final Parameter categoryParameter;
    private final PerBranchRewardMixtureCategoryDecoder categoryDecoder;
    private final RewardsMixtureBranchWeightProvider branchWeightProvider;

    public RewardMixtureCategoricalDiscontinuousPotentialProviderDynamic(
            final PerBranchRewardMixtureCategoryDecoder categoryDecoder,
            final RewardsAwareBranchModel rewardsAwareBranchModel,
            final TreeDataLikelihood treeDataLikelihood,
            final TreeDataLikelihood[] dependentTreeDataLikelihoods,
            final TreeDataLikelihood[] dependentContinuousTreeDataLikelihoods) {
        this(
                categoryDecoder,
                new RewardsMixtureBranchWeightProvider(
                        rewardsAwareBranchModel,
                        treeDataLikelihood,
                        dependentTreeDataLikelihoods,
                        dependentContinuousTreeDataLikelihoods));
    }

    public RewardMixtureCategoricalDiscontinuousPotentialProviderDynamic(
            final PerBranchRewardMixtureCategoryDecoder categoryDecoder,
            final RewardsMixtureBranchWeightProvider branchWeightProvider) {
        if (categoryDecoder == null) {
            throw new IllegalArgumentException("categoryDecoder must be non-null");
        }
        if (branchWeightProvider == null) {
            throw new IllegalArgumentException("branchWeightProvider must be non-null");
        }
        if (categoryDecoder.getCategoryParameter().getDimension() != branchWeightProvider.getBranchCount()) {
            throw new IllegalArgumentException(
                    "categoryParameter dimension must match branch count. Found " +
                            categoryDecoder.getCategoryParameter().getDimension() + " but expected " +
                            branchWeightProvider.getBranchCount());
        }

        this.categoryDecoder = categoryDecoder;
        this.categoryParameter = categoryDecoder.getCategoryParameter();
        this.branchWeightProvider = branchWeightProvider;
    }

    @Override
    public Parameter getParameter() {
        return categoryParameter;
    }

    @Override
    public int getDimension() {
        return categoryParameter.getDimension();
    }

    /**
     * Refreshes the decoder's per-branch insertion ranks once. Callers that
     * own the top of an HMC operator call (e.g.
     * MixedDiscontinuousHamiltonianMonteCarloOperator.doOperation) invoke
     * this exactly once before any coordinate integration; the per-crossing
     * methods below (getNextDiscontinuity, getPotentialDifferenceAcrossBoundary,
     * etc.) deliberately do NOT refresh, since they run many times per
     * operator call and refreshing there would reintroduce the O(branchCount)
     * per-crossing cost documented on PerBranchRewardMixtureCategoryDecoder.
     */
    @Override
    public void refresh() {
        categoryDecoder.refreshEmbedding();
    }

    @Override
    public double getLogDensity() {
        branchWeightProvider.refreshLikelihoodMessages();

        double logDensity = 0.0;
        for (int i = 0; i < getDimension(); i++) {
            final double contribution = getCurrentLogWeight(i);
            if (isNegativeInfinity(contribution)) {
                return Double.NEGATIVE_INFINITY;
            }
            logDensity += contribution;
        }
        return logDensity;
    }

    @Override
    public double getLogDensityAfterSingleCoordinateMove(final int index, final double proposedValue) {
        checkIndex(index);

        branchWeightProvider.refreshLikelihoodMessages();

        final double currentLogDensity = getLogDensityWithoutRefresh();
        if (isNegativeInfinity(currentLogDensity)) {
            return Double.NEGATIVE_INFINITY;
        }

        final RewardsMixtureBranchResamplingHelper.BranchWeights weights =
                branchWeightProvider.computeBranchWeightsForParameterIndex(index);
        final double currentContribution =
                logWeightForValue(weights, index, categoryParameter.getParameterValue(index));
        final double proposedContribution = logWeightForValue(weights, index, proposedValue);

        if (isNegativeInfinity(proposedContribution)) {
            return Double.NEGATIVE_INFINITY;
        }
        return currentLogDensity - currentContribution + proposedContribution;
    }

    @Override
    public double getPotentialDifference(final int index,
                                         final double currentValue,
                                         final double proposedValue) {
        checkIndex(index);

        final int currentCategory;
        final int proposedCategory;
        try {
            currentCategory = categoryDecoder.getCategoryForValue(index, currentValue);
        } catch (IllegalArgumentException e) {
            return Double.POSITIVE_INFINITY;
        }
        try {
            proposedCategory = categoryDecoder.getCategoryForValue(index, proposedValue);
        } catch (IllegalArgumentException e) {
            return Double.POSITIVE_INFINITY;
        }

        if (currentCategory == proposedCategory) {
            return 0.0;
        }

        branchWeightProvider.refreshLikelihoodMessages();
        final RewardsMixtureBranchResamplingHelper.BranchWeights weights =
                branchWeightProvider.computeBranchWeightsForParameterIndex(index);
        final double currentLogWeight =
                branchWeightProvider.getLogWeightForCategory(weights, currentCategory);
        final double proposedLogWeight =
                branchWeightProvider.getLogWeightForCategory(weights, proposedCategory);

        if (isNegativeInfinity(proposedLogWeight)) {
            return Double.POSITIVE_INFINITY;
        }
        if (isNegativeInfinity(currentLogWeight)) {
            return Double.NEGATIVE_INFINITY;
        }
        return currentLogWeight - proposedLogWeight;
    }

    /**
     * Mirrors RewardMixtureCategoricalDiscontinuousPotentialProvider's fast
     * override: only branch index's own before/after weight is computed,
     * not a full-tree sum -- see that class's javadoc for why (the O(branchCount)
     * redundant-recomputation bug fixed earlier this session).
     */
    @Override
    public double getPotentialDifferenceAcrossBoundary(final int index,
                                                        final double boundary,
                                                        final double direction) {
        checkIndex(index);
        if (direction == 0.0) {
            return 0.0;
        }

        final double before = Math.nextAfter(boundary, boundary - direction);
        final double after = Math.nextAfter(boundary, boundary + direction);

        final int currentCategory;
        final int proposedCategory;
        try {
            currentCategory = categoryDecoder.getCategoryForValue(index, before);
        } catch (IllegalArgumentException e) {
            return Double.NEGATIVE_INFINITY;
        }
        try {
            proposedCategory = categoryDecoder.getCategoryForValue(index, after);
        } catch (IllegalArgumentException e) {
            return Double.POSITIVE_INFINITY;
        }

        if (currentCategory == proposedCategory) {
            return 0.0;
        }

        branchWeightProvider.refreshLikelihoodMessages();
        final RewardsMixtureBranchResamplingHelper.BranchWeights weights =
                branchWeightProvider.computeBranchWeightsForParameterIndex(index);
        final double currentLogWeight =
                branchWeightProvider.getLogWeightForCategory(weights, currentCategory);
        final double proposedLogWeight =
                branchWeightProvider.getLogWeightForCategory(weights, proposedCategory);

        if (isNegativeInfinity(proposedLogWeight)) {
            return Double.POSITIVE_INFINITY;
        }
        if (isNegativeInfinity(currentLogWeight)) {
            return Double.NEGATIVE_INFINITY;
        }
        return currentLogWeight - proposedLogWeight;
    }

    @Override
    public double getNextDiscontinuity(final int index,
                                       final double currentValue,
                                       final double direction) {
        checkIndex(index);
        return categoryDecoder.getNextBoundary(currentValue, direction);
    }

    public int getRewardMixtureCategory(final int index) {
        checkIndex(index);
        categoryDecoder.refreshEmbedding();
        return categoryDecoder.getCategoryForParameterIndex(index);
    }

    public int getCategoryCount() {
        return categoryDecoder.getCategoryCount();
    }

    private double getLogDensityWithoutRefresh() {
        double logDensity = 0.0;
        for (int i = 0; i < getDimension(); i++) {
            final double contribution = getCurrentLogWeight(i);
            if (isNegativeInfinity(contribution)) {
                return Double.NEGATIVE_INFINITY;
            }
            logDensity += contribution;
        }
        return logDensity;
    }

    private double getCurrentLogWeight(final int index) {
        final RewardsMixtureBranchResamplingHelper.BranchWeights weights =
                branchWeightProvider.computeBranchWeightsForParameterIndex(index);
        return logWeightForValue(weights, index, categoryParameter.getParameterValue(index));
    }

    private double logWeightForValue(final RewardsMixtureBranchResamplingHelper.BranchWeights weights,
                                     final int index,
                                     final double value) {
        try {
            return branchWeightProvider.getLogWeightForCategory(
                    weights,
                    categoryDecoder.getCategoryForValue(index, value));
        } catch (IllegalArgumentException e) {
            return Double.NEGATIVE_INFINITY;
        }
    }

    private void checkIndex(final int index) {
        if (index < 0 || index >= getDimension()) {
            throw new IllegalArgumentException("Coordinate index out of bounds: " + index);
        }
    }

    private static boolean isNegativeInfinity(final double value) {
        return Double.isInfinite(value) && value < 0.0;
    }
}
