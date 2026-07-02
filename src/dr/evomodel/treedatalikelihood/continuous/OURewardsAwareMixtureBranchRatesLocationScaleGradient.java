package dr.evomodel.treedatalikelihood.continuous;

import dr.evolution.tree.NodeRef;
import dr.evomodel.branchratemodel.RewardsAwareMixtureBranchRates;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.preorder.BranchSufficientStatistics;
import dr.inference.model.Parameter;

public final class OURewardsAwareMixtureBranchRatesLocationScaleGradient
        extends AbstractRewardsAwareMixtureBranchRatesLocationScaleGradient {

    private final OUBranchTimeGradient branchTimeGradient;

    public OURewardsAwareMixtureBranchRatesLocationScaleGradient(
            String traitName,
            TreeDataLikelihood treeDataLikelihood,
            RewardsAwareMixtureBranchRates branchRateModel,
            Parameter locationScaleParameter,
            Double tolerance) {

        super("OU dependent location-scale gradient",
                traitName,
                treeDataLikelihood,
                branchRateModel,
                locationScaleParameter,
                tolerance);

        this.branchTimeGradient = new OUBranchTimeGradient(
                likelihoodDelegate.getTraitDim(),
                tree,
                likelihoodDelegate);
    }

    @Override
    protected void accumulateBranchGradient(
            double[] result,
            NodeRef node,
            BranchSufficientStatistics statistics,
            double rawReward,
            double rate) {

        final double scaledTime = getScaledTime(node);
        final double branchGradient = branchTimeGradient.getGradientForBranch(statistics, node)[0];

        result[LOCATION_INDEX] += branchGradient * scaledTime / rate;
        result[SCALE_INDEX] += branchGradient * scaledTime * rawReward / rate;
    }

    private double getScaledTime(NodeRef node) {
        final int matrixIndex = likelihoodDelegate.getDiffusionProcessDelegate()
                .getMatrixBufferOffsetIndex(node.getNumber());
        return likelihoodDelegate.getIntegrator().getBranchLength(matrixIndex);
    }
}
