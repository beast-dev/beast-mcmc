package dr.inference.operators;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.RewardsAwareMixtureBranchRates;
import dr.evomodel.treedatalikelihood.DataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.inference.model.Parameter;

/**
 * Exact continuous-trait dependent evidence for candidate reward values.
 *
 * This intentionally uses the continuous likelihood itself as the oracle. It is
 * slower than a branch-local Gaussian contraction, but supports the existing BM
 * and OU likelihood delegates immediately and gives an exact reference for a
 * later local evaluator.
 */
public final class ContinuousRewardDependentEdgeEvidenceProvider implements RewardDependentEdgeEvidenceProvider {

    private final TreeDataLikelihood treeDataLikelihood;
    private final Tree tree;
    private final RewardsAwareMixtureBranchRates rewardBranchRates;
    private final Parameter indicator;
    private final Parameter continuousRewards;

    public ContinuousRewardDependentEdgeEvidenceProvider(final TreeDataLikelihood treeDataLikelihood) {
        if (treeDataLikelihood == null) {
            throw new IllegalArgumentException("treeDataLikelihood must be non-null");
        }

        final DataLikelihoodDelegate delegate = treeDataLikelihood.getDataLikelihoodDelegate();
        if (!(delegate instanceof ContinuousDataLikelihoodDelegate)) {
            throw new IllegalArgumentException(
                    "Continuous reward evidence requires ContinuousDataLikelihoodDelegate, found " +
                            (delegate == null ? "null" : delegate.getClass().getName())
            );
        }

        final BranchRateModel branchRateModel = treeDataLikelihood.getBranchRateModel();
        if (!(branchRateModel instanceof RewardsAwareMixtureBranchRates)) {
            throw new IllegalArgumentException(
                    "Continuous reward evidence requires RewardsAwareMixtureBranchRates, found " +
                            branchRateModel.getClass().getName()
            );
        }

        this.treeDataLikelihood = treeDataLikelihood;
        this.tree = treeDataLikelihood.getTree();
        this.rewardBranchRates = (RewardsAwareMixtureBranchRates) branchRateModel;
        this.indicator = rewardBranchRates.getIndicator();
        this.continuousRewards = rewardBranchRates.getRateParameter();
    }

    @Override
    public void prepare() {
        treeDataLikelihood.makeDirty();
    }

    @Override
    public double logEvidence(final int branchNodeNumber, final double rawReward) {
        if (branchNodeNumber < 0 || branchNodeNumber >= tree.getNodeCount()) {
            throw new IllegalArgumentException("branchNodeNumber out of range: " + branchNodeNumber);
        }
        final NodeRef node = tree.getNode(branchNodeNumber);
        if (tree.isRoot(node)) {
            throw new IllegalArgumentException("Root node has no branch: " + branchNodeNumber);
        }

        final int parameterIndex = rewardBranchRates.getParameterIndexFromNode(node);
        if (parameterIndex < 0 || parameterIndex >= continuousRewards.getDimension()) {
            throw new IllegalArgumentException(
                    "Invalid continuous reward parameter index " + parameterIndex +
                            " for branch node " + branchNodeNumber);
        }

        final double oldReward = continuousRewards.getParameterValue(parameterIndex);
        final double oldIndicator = indicator.getParameterValue(parameterIndex);

        try {
            continuousRewards.setParameterValueQuietly(parameterIndex, rawReward);
            indicator.setParameterValueQuietly(parameterIndex, 0.0);
            treeDataLikelihood.makeDirty();
            final double logLikelihood = treeDataLikelihood.getLogLikelihood();
            return Double.isFinite(logLikelihood) ? logLikelihood : Double.NEGATIVE_INFINITY;
        } finally {
            continuousRewards.setParameterValueQuietly(parameterIndex, oldReward);
            indicator.setParameterValueQuietly(parameterIndex, oldIndicator);
            treeDataLikelihood.makeDirty();
        }
    }
}
