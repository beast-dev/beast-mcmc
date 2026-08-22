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
 * Deliberately excludes value-only decoding (decode a raw axis value with no
 * branch context) and cut-boundary lookups: both are branch-dependent for
 * the per-branch decoder, so callers that need them (currently only
 * {@code RewardMixtureCategoricalDiscontinuousPotentialProvider}'s
 * boundary-crossing evaluation) use the concrete decoder type directly
 * rather than forcing every implementation to accept a branch index it may
 * not need. {@link #getNextBoundary} is the one value-based method safe to
 * share: the axis itself (cut values) is identical for every branch, so
 * finding the next boundary from a raw value never needs a branch index.
 */
public interface RewardMixtureCategoryDecoding {

    Parameter getCategoryParameter();

    Parameter getCutParameter();

    int getCategoryCount();

    void refreshEmbedding();

    int getCategoryForParameterIndex(int parameterIndex);

    double getNextBoundary(double value, double direction);

    boolean isAtomic(int parameterIndex);

    int getAtomicState(int parameterIndex);
}
