package dr.evomodel.branchratemodel;

import dr.inference.hmc.EmbeddedOrdinalParameter;
import dr.inference.model.Parameter;

/**
 * Branch-specific reward-mixture category decoder: atomic states are
 * ordered along the embedded axis by reward value, and the continuous
 * category sits between whichever two atomic states currently bracket that
 * branch's live total.rewards.cts value (see
 * REWARD_CATEGORY_DYNAMIC_ORDERING_PLAN.md).
 *
 * The axis itself (cut values 0,1,...,K+1) is identical for every branch and
 * built once, exactly like {@link RewardMixtureCategoryDecoder}'s shared
 * embedding -- only the bucket-to-category translation is branch-specific,
 * via a per-branch insertion rank. That translation is computed for every
 * branch once per {@link #refreshEmbedding()} call (O(branchCount log K))
 * and then reused for every decode call until the next refresh, so decode
 * itself stays O(log K) -- the same asymptotic cost as the shared decoder,
 * not O(branchCount). Callers MUST only invoke refreshEmbedding() once per
 * HMC operator call (matching where the shared decoder's refresh already
 * happens today), never per boundary crossing within one trajectory: that is
 * both a correctness requirement (the discontinuous integrator needs a
 * static potential landscape for the duration of one integrated trajectory)
 * and what keeps this decoder's cost light rather than repeating
 * O(branchCount) work per crossing -- the same class of bug fixed in
 * RewardMixtureCategoricalDiscontinuousPotentialProvider earlier this
 * session.
 */
public final class PerBranchRewardMixtureCategoryDecoder implements RewardMixtureCategoryDecoding {

    private final Parameter categoryParameter;
    private final Parameter cutParameter;
    private final Parameter ctsParameter;
    private final RewardRates rewardRates;
    private final int atomicStateCount;

    private EmbeddedOrdinalParameter embedding;
    private PerBranchCategoryLayout layout;
    private int[] insertionRankByParameterIndex;

    public PerBranchRewardMixtureCategoryDecoder(final Parameter categoryParameter,
                                                 final Parameter cutParameter,
                                                 final Parameter ctsParameter,
                                                 final RewardRates rewardRates,
                                                 final int atomicStateCount,
                                                 final int expectedDimension) {
        if (categoryParameter == null) {
            throw new IllegalArgumentException("categoryParameter must be non-null");
        }
        if (cutParameter == null) {
            throw new IllegalArgumentException("cutParameter must be non-null");
        }
        if (ctsParameter == null) {
            throw new IllegalArgumentException("ctsParameter must be non-null");
        }
        if (rewardRates == null) {
            throw new IllegalArgumentException("rewardRates must be non-null");
        }
        if (atomicStateCount < 1) {
            throw new IllegalArgumentException("atomicStateCount must be positive");
        }
        if (categoryParameter.getDimension() != expectedDimension) {
            throw new IllegalArgumentException(
                    "categoryParameter dimension must be " + expectedDimension +
                            " but is " + categoryParameter.getDimension());
        }
        if (ctsParameter.getDimension() != expectedDimension) {
            throw new IllegalArgumentException(
                    "ctsParameter dimension must be " + expectedDimension +
                            " but is " + ctsParameter.getDimension());
        }
        if (cutParameter.getDimension() != atomicStateCount + 2) {
            throw new IllegalArgumentException(
                    "cutParameter dimension must be atomicStateCount + 2 (" +
                            (atomicStateCount + 2) + ") but is " + cutParameter.getDimension());
        }

        this.categoryParameter = categoryParameter;
        this.cutParameter = cutParameter;
        this.ctsParameter = ctsParameter;
        this.rewardRates = rewardRates;
        this.atomicStateCount = atomicStateCount;
        refreshEmbedding();
    }

    @Override
    public Parameter getCategoryParameter() {
        return categoryParameter;
    }

    @Override
    public Parameter getCutParameter() {
        return cutParameter;
    }

    @Override
    public int getCategoryCount() {
        return atomicStateCount + 1;
    }

    @Override
    public void refreshEmbedding() {
        embedding = new EmbeddedOrdinalParameter(cutParameter.getParameterValues());
        if (embedding.getStateCount() != getCategoryCount()) {
            throw new IllegalArgumentException(
                    "Embedding state count must be " + getCategoryCount() +
                            " but is " + embedding.getStateCount());
        }

        final double[] rewardByState = new double[atomicStateCount];
        for (int s = 0; s < atomicStateCount; s++) {
            rewardByState[s] = rewardRates.getRawReward(s);
        }
        layout = PerBranchCategoryLayout.fromAtomicRewards(rewardByState);

        final int dim = ctsParameter.getDimension();
        if (insertionRankByParameterIndex == null || insertionRankByParameterIndex.length != dim) {
            insertionRankByParameterIndex = new int[dim];
        }
        for (int p = 0; p < dim; p++) {
            insertionRankByParameterIndex[p] = layout.insertionRank(ctsParameter.getParameterValue(p));
        }
    }

    @Override
    public int getCategoryForParameterIndex(final int parameterIndex) {
        checkParameterIndex(parameterIndex);
        final int bucket = embedding.getStateIndex(categoryParameter.getParameterValue(parameterIndex));
        return layout.categoryForBucket(insertionRankByParameterIndex[parameterIndex], bucket);
    }

    /** Decode a raw embedded-coordinate value for a specific branch (unlike the shared decoder, which axis bucket is "continuous" depends on that branch's insertion rank). */
    public int getCategoryForValue(final int parameterIndex, final double value) {
        checkParameterIndex(parameterIndex);
        final int bucket = embedding.getStateIndex(value);
        return layout.categoryForBucket(insertionRankByParameterIndex[parameterIndex], bucket);
    }

    @Override
    public double getNextBoundary(final double value, final double direction) {
        // The axis (cut values) is identical for every branch, so finding the
        // next boundary from a raw value never needs branch context.
        return embedding.getNextBoundary(value, direction);
    }

    @Override
    public boolean isAtomic(final int parameterIndex) {
        return RewardMixtureCategoryDecoder.isAtomicCategory(getCategoryForParameterIndex(parameterIndex));
    }

    @Override
    public int getAtomicState(final int parameterIndex) {
        return RewardMixtureCategoryDecoder.getAtomicStateForCategory(getCategoryForParameterIndex(parameterIndex));
    }

    @Override
    public double getLowerCut(final int parameterIndex, final int category) {
        checkParameterIndex(parameterIndex);
        return embedding.getLowerCut(layout.bucketForCategory(insertionRankByParameterIndex[parameterIndex], category));
    }

    @Override
    public double getUpperCut(final int parameterIndex, final int category) {
        checkParameterIndex(parameterIndex);
        return embedding.getUpperCut(layout.bucketForCategory(insertionRankByParameterIndex[parameterIndex], category));
    }

    public double getLowerCutForCategoryAtCtsValue(final int parameterIndex,
                                                   final int category,
                                                   final double ctsValue) {
        checkParameterIndex(parameterIndex);
        final int rank = layout.insertionRank(ctsValue);
        return embedding.getLowerCut(layout.bucketForCategory(rank, category));
    }

    public double getUpperCutForCategoryAtCtsValue(final int parameterIndex,
                                                   final int category,
                                                   final double ctsValue) {
        checkParameterIndex(parameterIndex);
        final int rank = layout.insertionRank(ctsValue);
        return embedding.getUpperCut(layout.bucketForCategory(rank, category));
    }

    private void checkParameterIndex(final int parameterIndex) {
        if (parameterIndex < 0 || parameterIndex >= categoryParameter.getDimension()) {
            throw new IllegalArgumentException("Category parameter index out of range: " + parameterIndex);
        }
    }
}
