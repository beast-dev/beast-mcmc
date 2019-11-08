package dr.evomodel.treedatalikelihood.discrete;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.BranchSpecificFixedEffects;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;

/**
 * @author Marc A. Suchard
 * @author Xiang Ji
 */
public class LocationGradient extends HyperParameterBranchRateGradient {

    private final BranchSpecificFixedEffects fixedEffects;

    public LocationGradient(String traitName, TreeDataLikelihood treeDataLikelihood,
                            BeagleDataLikelihoodDelegate likelihoodDelegate,
                            BranchSpecificFixedEffects fixedEffects,
                            boolean useHessian) {

        super(traitName, treeDataLikelihood, likelihoodDelegate, fixedEffects.getFixedEffectsParameter(), useHessian);

        this.fixedEffects = fixedEffects;
    }

    @Override
    double[] getDifferential(Tree tree, NodeRef node) {

        double rate = branchRateModel.getBranchRate(tree, node);
        double[] results = fixedEffects.getDifferential(rate, tree, node);

        return results;
    }
}
