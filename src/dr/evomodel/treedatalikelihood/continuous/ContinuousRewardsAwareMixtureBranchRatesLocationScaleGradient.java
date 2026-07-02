package dr.evomodel.treedatalikelihood.continuous;

import dr.evolution.tree.NodeRef;
import dr.evomodel.branchratemodel.RewardsAwareMixtureBranchRates;
import dr.evomodel.continuous.SparseBandedMultivariateDiffusionModel;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.preorder.BranchSufficientStatistics;
import dr.inference.model.Parameter;

/**
 * Continuous-trait gradient with respect to the global location and scale in a
 * linear location-scale RewardsAwareMixtureBranchRates transform.
 */
public final class ContinuousRewardsAwareMixtureBranchRatesLocationScaleGradient
        extends AbstractRewardsAwareMixtureBranchRatesLocationScaleGradient {

    private final BranchRateGradient.ContinuousTraitGradientForBranch branchProvider;

    public ContinuousRewardsAwareMixtureBranchRatesLocationScaleGradient(
            String traitName,
            TreeDataLikelihood treeDataLikelihood,
            RewardsAwareMixtureBranchRates branchRateModel,
            Parameter locationScaleParameter,
            Double tolerance) {

        super("Continuous dependent location-scale gradient",
                traitName,
                treeDataLikelihood,
                branchRateModel,
                locationScaleParameter,
                tolerance);

        final int dim = likelihoodDelegate.getTraitDim();
        if (likelihoodDelegate.getDiffusionModel() instanceof SparseBandedMultivariateDiffusionModel) {
            branchProvider = new BranchRateGradient.ContinuousTraitGradientForBranch.Sparse(dim);
        } else {
            branchProvider = new BranchRateGradient.ContinuousTraitGradientForBranch.Dense(dim);
        }
    }

    @Override
    protected void accumulateBranchGradient(
            double[] result,
            NodeRef node,
            BranchSufficientStatistics statistics,
            double rawReward,
            double rate) {

        result[LOCATION_INDEX] += branchProvider.getGradientForBranch(statistics, 1.0 / rate);
        result[SCALE_INDEX] += branchProvider.getGradientForBranch(statistics, rawReward / rate);
    }
}
