package dr.evomodel.branchratemodel;

import dr.inference.model.Parameter;

/**
 * RewardMixtureBranchRateModel specialization for the embedded-categorical
 * state representation (as opposed to the legacy indicator/atom-index
 * representation), implemented by both {@link RewardsAwareCategoricalMixtureBranchRates}
 * and {@link RewardsAwareCategoricalMixtureBranchRatesDynamic}. Lets code that
 * only needs the category parameter/cuts (not the branch-rate dispatch logic
 * itself) work with either implementation.
 */
public interface RewardMixtureCategoricalBranchRateModel extends RewardMixtureBranchRateModel {

    Parameter getCategoryParameter();

    Parameter getCategoryCutParameter();
}
