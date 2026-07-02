package dr.evomodel.treedatalikelihood.continuous;

import dr.evolution.tree.NodeRef;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.preorder.BranchSufficientStatistics;
import dr.inference.model.Parameter;

public final class OUBranchRateGradient extends AbstractContinuousBranchRateGradient {

    private final OUBranchTimeGradient branchTimeGradient;

    public OUBranchRateGradient(String traitName,
                                TreeDataLikelihood treeDataLikelihood,
                                ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                Parameter rateParameter) {

        super("OU branch-rate gradient (branch-time derivative)",
                traitName,
                treeDataLikelihood,
                likelihoodDelegate,
                rateParameter);

        this.branchTimeGradient = new OUBranchTimeGradient(
                likelihoodDelegate.getTraitDim(),
                tree,
                likelihoodDelegate);
    }

    @Override
    protected double getGradientWrtBranchTime(BranchSufficientStatistics statistics, NodeRef node) {
        return branchTimeGradient.getGradientForBranch(statistics, node)[0];
    }
}
