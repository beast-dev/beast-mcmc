package dr.inference.hmc;

import dr.evomodel.branchmodel.RewardsAwareBranchModel;
import dr.evomodel.branchmodel.RewardMixtureAtomicPseudoPrior;
import dr.evomodel.branchratemodel.RewardMixtureCategoryDecoder;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.inference.model.Parameter;
import dr.inference.operators.RewardMixturePerformanceStats;
import dr.inference.operators.RewardsMixtureBranchResamplingHelper;
import dr.inference.operators.RewardsMixtureBranchWeightProvider;

/**
 * Discontinuous-HMC potential provider for reward-mixture branch states.
 *
 * A single embedded categorical coordinate represents each branch:
 * category 0 is the continuous reward state, and categories 1..K are atomic
 * no-jump states 0..K-1.
 */
public final class RewardMixtureCategoricalDiscontinuousPotentialProvider
        implements DiscontinuousPotentialProvider {

    private final Parameter categoryParameter;
    private final RewardMixtureCategoryDecoder categoryDecoder;
    private final RewardsMixtureBranchWeightProvider branchWeightProvider;
    private final RewardMixtureAtomicPseudoPrior atomicPseudoPrior;

    public RewardMixtureCategoricalDiscontinuousPotentialProvider(
            final Parameter categoryParameter,
            final Parameter categoryCuts,
            final RewardsAwareBranchModel rewardsAwareBranchModel,
            final TreeDataLikelihood treeDataLikelihood,
            final TreeDataLikelihood[] dependentTreeDataLikelihoods,
            final TreeDataLikelihood[] dependentContinuousTreeDataLikelihoods) {
        this(categoryParameter,
                categoryCuts,
                rewardsAwareBranchModel,
                treeDataLikelihood,
                dependentTreeDataLikelihoods,
                dependentContinuousTreeDataLikelihoods,
                null);
    }

    public RewardMixtureCategoricalDiscontinuousPotentialProvider(
            final Parameter categoryParameter,
            final Parameter categoryCuts,
            final RewardsAwareBranchModel rewardsAwareBranchModel,
            final TreeDataLikelihood treeDataLikelihood,
            final TreeDataLikelihood[] dependentTreeDataLikelihoods,
            final TreeDataLikelihood[] dependentContinuousTreeDataLikelihoods,
            final RewardMixtureAtomicPseudoPrior atomicPseudoPrior) {
        this(
                categoryParameter,
                categoryCuts,
                new RewardsMixtureBranchWeightProvider(
                        rewardsAwareBranchModel,
                        treeDataLikelihood,
                        dependentTreeDataLikelihoods,
                        dependentContinuousTreeDataLikelihoods),
                atomicPseudoPrior);
    }

    public RewardMixtureCategoricalDiscontinuousPotentialProvider(
            final Parameter categoryParameter,
            final Parameter categoryCuts,
            final RewardsMixtureBranchWeightProvider branchWeightProvider) {
        this(categoryParameter, categoryCuts, branchWeightProvider, null);
    }

    public RewardMixtureCategoricalDiscontinuousPotentialProvider(
            final Parameter categoryParameter,
            final Parameter categoryCuts,
            final RewardsMixtureBranchWeightProvider branchWeightProvider,
            final RewardMixtureAtomicPseudoPrior atomicPseudoPrior) {
        if (categoryParameter == null) {
            throw new IllegalArgumentException("categoryParameter must be non-null");
        }
        if (branchWeightProvider == null) {
            throw new IllegalArgumentException("branchWeightProvider must be non-null");
        }
        if (categoryParameter.getDimension() != branchWeightProvider.getBranchCount()) {
            throw new IllegalArgumentException(
                    "categoryParameter dimension must match branch count. Found " +
                            categoryParameter.getDimension() + " but expected " +
                            branchWeightProvider.getBranchCount());
        }

        this.categoryParameter = categoryParameter;
        this.branchWeightProvider = branchWeightProvider;
        this.atomicPseudoPrior = atomicPseudoPrior;
        this.categoryDecoder = new RewardMixtureCategoryDecoder(
                categoryParameter,
                categoryCuts,
                branchWeightProvider.getStateCount(),
                branchWeightProvider.getBranchCount());
    }

    @Override
    public Parameter getParameter() {
        return categoryParameter;
    }

    @Override
    public int getDimension() {
        return categoryParameter.getDimension();
    }

    @Override
    public void refresh() {
        categoryDecoder.refreshEmbedding();
        branchWeightProvider.beginOperationCache();
    }

    @Override
    public void refreshAfterPositionUpdate() {
        categoryDecoder.refreshEmbedding();
        branchWeightProvider.clearOperationCache(
                RewardMixturePerformanceStats.OperationCacheClearReason.FINAL_POSITION_REFRESH);
    }

    @Override
    public void clearOperationCache() {
        clearOperationCache(RewardMixturePerformanceStats.OperationCacheClearReason.UNKNOWN);
    }

    @Override
    public void clearOperationCache(final RewardMixturePerformanceStats.OperationCacheClearReason reason) {
        branchWeightProvider.clearOperationCache(reason);
    }

    @Override
    public double getLogDensity() {
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

        final double currentLogDensity = getLogDensityWithoutRefresh();
        if (isNegativeInfinity(currentLogDensity)) {
            return Double.NEGATIVE_INFINITY;
        }

        final RewardsMixtureBranchResamplingHelper.BranchWeights weights = getBranchWeights(index);
        final double currentContribution =
                logWeightForValue(index, weights, categoryParameter.getParameterValue(index));
        final double proposedContribution = logWeightForValue(index, weights, proposedValue);

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
            currentCategory = categoryDecoder.getCategoryForValue(currentValue);
        } catch (IllegalArgumentException e) {
            return Double.POSITIVE_INFINITY;
        }
        try {
            proposedCategory = categoryDecoder.getCategoryForValue(proposedValue);
        } catch (IllegalArgumentException e) {
            return Double.POSITIVE_INFINITY;
        }

        if (currentCategory == proposedCategory) {
            return 0.0;
        }

        final RewardsMixtureBranchResamplingHelper.BranchWeights weights = getBranchWeights(index);
        final double currentLogWeight =
                getLogWeightForCategory(index, weights, currentCategory);
        final double proposedLogWeight =
                getLogWeightForCategory(index, weights, proposedCategory);

        if (isNegativeInfinity(proposedLogWeight)) {
            return Double.POSITIVE_INFINITY;
        }
        if (isNegativeInfinity(currentLogWeight)) {
            return Double.NEGATIVE_INFINITY;
        }
        return currentLogWeight - proposedLogWeight;
    }

    /**
     * Overrides the default interface implementation, which evaluates two full
     * {@link #getLogDensityAfterSingleCoordinateMove} calls -- each summing every
     * branch's contribution via {@code getLogDensityWithoutRefresh} -- even though
     * only branch {@code index}'s own weight differs across a boundary crossing and
     * every other branch's contribution cancels in the caller's subtraction. This
     * mirrors {@link #getPotentialDifference}: touch only branch {@code index}.
     * At small branch counts the O(branchCount) waste in the default path is
     * invisible; at real tree sizes (hundreds of branches) it turns every single
     * boundary crossing the discontinuous coordinate integrator evaluates into an
     * O(branchCount) pass, compounding across all crossings in a step.
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
            currentCategory = categoryDecoder.getCategoryForValue(before);
        } catch (IllegalArgumentException e) {
            return Double.NEGATIVE_INFINITY;
        }
        try {
            proposedCategory = categoryDecoder.getCategoryForValue(after);
        } catch (IllegalArgumentException e) {
            return Double.POSITIVE_INFINITY;
        }

        if (currentCategory == proposedCategory) {
            return 0.0;
        }

        final RewardsMixtureBranchResamplingHelper.BranchWeights weights = getBranchWeights(index);
        final double currentLogWeight =
                getLogWeightForCategory(index, weights, currentCategory);
        final double proposedLogWeight =
                getLogWeightForCategory(index, weights, proposedCategory);

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
        final RewardsMixtureBranchResamplingHelper.BranchWeights weights = getBranchWeights(index);
        return logWeightForValue(index, weights, categoryParameter.getParameterValue(index));
    }

    private RewardsMixtureBranchResamplingHelper.BranchWeights getBranchWeights(final int index) {
        return branchWeightProvider.getOperationCachedBranchWeightsForParameterIndex(index);
    }

    private double logWeightForValue(final int index,
                                     final RewardsMixtureBranchResamplingHelper.BranchWeights weights,
                                     final double value) {
        try {
            return getLogWeightForCategory(
                    index,
                    weights,
                    categoryDecoder.getCategoryForValue(value));
        } catch (IllegalArgumentException e) {
            return Double.NEGATIVE_INFINITY;
        }
    }

    private double getLogWeightForCategory(final int index,
                                           final RewardsMixtureBranchResamplingHelper.BranchWeights weights,
                                           final int category) {
        return branchWeightProvider.getLogWeightForCategory(weights, category) +
                getPseudoPriorLogDensity(index, category);
    }

    private double getPseudoPriorLogDensity(final int index, final int category) {
        return atomicPseudoPrior == null ? 0.0 : atomicPseudoPrior.getLogDensityForCategory(index, category);
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
