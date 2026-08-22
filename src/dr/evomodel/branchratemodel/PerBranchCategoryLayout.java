package dr.evomodel.branchratemodel;

import java.util.Arrays;

/**
 * Pure logic for laying out reward-mixture categories along an embedded
 * coordinate axis so atomic states are ordered by reward value and the
 * continuous category sits between whichever two atomic states currently
 * bracket a branch's live total.rewards.cts value.
 *
 * "Bucket" is the raw 0..K axis position produced by EmbeddedOrdinalParameter
 * (whose cut values 0,1,...,K+1 never change and are shared by every
 * branch); "category" is the physical identity (0 = continuous, k = atomic
 * state k-1) that the rest of the reward-mixture code already understands.
 * Only the bucket-to-category translation is branch-specific (via that
 * branch's insertion rank); everything about the axis itself -- cut values,
 * next-boundary search -- stays shared and untouched.
 *
 * No Parameter/Model dependencies: this class is deliberately pure so its
 * correctness (including the tie/edge cases) can be tested exhaustively
 * without any BEAST scaffolding.
 */
public final class PerBranchCategoryLayout {

    public static final int CONTINUOUS_CATEGORY = 0;

    private final double[] sortedRewards;      // ascending, length K
    private final int[] sortedStateIndices;    // atomic state at each sorted rank, length K
    private final int[] rankOfState;           // inverse of sortedStateIndices, length K

    private PerBranchCategoryLayout(final double[] sortedRewards,
                                    final int[] sortedStateIndices,
                                    final int[] rankOfState) {
        this.sortedRewards = sortedRewards;
        this.sortedStateIndices = sortedStateIndices;
        this.rankOfState = rankOfState;
    }

    /**
     * @param rewardByState raw reward value for atomic state i at rewardByState[i]
     */
    public static PerBranchCategoryLayout fromAtomicRewards(final double[] rewardByState) {
        final int k = rewardByState.length;
        if (k < 1) {
            throw new IllegalArgumentException("At least one atomic state is required");
        }

        final Integer[] order = new Integer[k];
        for (int i = 0; i < k; i++) {
            order[i] = i;
        }
        Arrays.sort(order, (a, b) -> Double.compare(rewardByState[a], rewardByState[b]));

        final double[] sortedRewards = new double[k];
        final int[] sortedStateIndices = new int[k];
        final int[] rankOfState = new int[k];
        for (int rank = 0; rank < k; rank++) {
            final int state = order[rank];
            sortedRewards[rank] = rewardByState[state];
            sortedStateIndices[rank] = state;
            rankOfState[state] = rank;
        }
        return new PerBranchCategoryLayout(sortedRewards, sortedStateIndices, rankOfState);
    }

    public int getAtomicStateCount() {
        return sortedRewards.length;
    }

    public double getSortedReward(final int rank) {
        return sortedRewards[rank];
    }

    public int getSortedStateIndex(final int rank) {
        return sortedStateIndices[rank];
    }

    /**
     * Number of atomic rewards strictly less than ctsValue -- the axis bucket
     * the continuous category should occupy for a branch whose current
     * total.rewards.cts value is ctsValue. O(log K). A ctsValue exactly
     * equal to an atomic reward counts as "not less than" that reward, so
     * ties place the continuous slot on the lower side deterministically
     * (see class-level note on the 2026-08-22 total.rewards.cts=0.5
     * coincidence this mirrors but does not itself resolve).
     */
    public int insertionRank(final double ctsValue) {
        int lo = 0;
        int hi = sortedRewards.length;
        while (lo < hi) {
            final int mid = (lo + hi) >>> 1;
            if (sortedRewards[mid] < ctsValue) {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }
        return lo;
    }

    /** Axis bucket (0..K) that the given category occupies for a branch with this insertion rank. */
    public int bucketForCategory(final int rank, final int category) {
        if (category == CONTINUOUS_CATEGORY) {
            return rank;
        }
        final int atomicState = category - 1;
        final int stateRank = rankOfState[atomicState];
        return stateRank < rank ? stateRank : stateRank + 1;
    }

    /** Category occupying the given axis bucket for a branch with this insertion rank. */
    public int categoryForBucket(final int rank, final int bucket) {
        if (bucket == rank) {
            return CONTINUOUS_CATEGORY;
        }
        final int stateRank = bucket < rank ? bucket : bucket - 1;
        return sortedStateIndices[stateRank] + 1;
    }
}
