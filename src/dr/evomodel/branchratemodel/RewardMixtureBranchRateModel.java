package dr.evomodel.branchratemodel;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.inference.model.Parameter;

/**
 * Shared branch-rate surface for reward-mixture models.
 *
 * The selected branch rate may be continuous or atomic, but resampling and
 * discontinuous-HMC proposals also need access to the latent continuous reward
 * associated with a branch.
 */
public interface RewardMixtureBranchRateModel {

    Parameter getRateParameter();

    int getParameterIndexFromNode(NodeRef node);

    double getBranchRate(Tree tree, NodeRef node);

    double getBranchRateForRawReward(Tree tree, NodeRef node, double rawReward);

    double getContinuousRawReward(Tree tree, NodeRef node);

    double getRawRewardForAtomState(int stateIndex);

    RewardRates getRewardRates();
}
