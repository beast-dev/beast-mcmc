package dr.evomodel.branchratemodel;

import dr.inference.model.Parameter;

/**
 * Common surface shared by {@link RewardMixtureCategoryDecoder} (one shared
 * embedding for the whole model) and {@code PerBranchRewardMixtureCategoryDecoder}
 * (branch-specific embedding, continuous slot tracks that branch's live
 * total.rewards.cts value). Consumers that only need "which category is this
 * branch in," keyed by parameter index, can depend on this interface and
 * work unchanged with either decoder.
 *
 * Cut-boundary lookups require a branch index because the physical category
 * occupying an embedded-axis bucket is branch-dependent for the per-branch
 * decoder. Static decoders ignore the branch index after bounds checking.
 * {@link #getNextBoundary} is still value-only because the axis itself (cut
 * values) is identical for every branch.
 */
public interface RewardMixtureCategoryDecoding {

    Parameter getCategoryParameter();

    Parameter getCutParameter();

    int getCategoryCount();

    void refreshEmbedding();

    int getCategoryForParameterIndex(int parameterIndex);

    double getLowerCut(int parameterIndex, int category);

    double getUpperCut(int parameterIndex, int category);

    double getNextBoundary(double value, double direction);

    boolean isAtomic(int parameterIndex);

    int getAtomicState(int parameterIndex);
}
