package dr.evomodel.treedatalikelihood.discrete;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.BranchSpecificFixedEffects;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;

/**
 * @author Marc A. Suchard
 */
public class NewLocationGradient extends HyperParameterBranchRateGradient {

    private final BranchSpecificFixedEffects fixedEffects;

    public NewLocationGradient(String traitName, TreeDataLikelihood treeDataLikelihood,
                               BeagleDataLikelihoodDelegate likelihoodDelegate,
                               BranchSpecificFixedEffects fixedEffects,
                               boolean useHessian) {

        super(traitName, treeDataLikelihood, likelihoodDelegate, fixedEffects.getFixedEffectsParameter(), useHessian);

        this.fixedEffects = fixedEffects;
    }

    @Override
    double[] getDifferential(Tree tree, NodeRef node) {
        double rate = branchRateModel.getBranchRate(tree, node);

        double[] result = fixedEffects.getDesignVector(tree, node);

        final double multiplier = rate / fixedEffects.getEffect(tree, node);

        for (int i = 0; i < result.length; i++) {
            result[i] *= multiplier;
        }

        return result;
    }
}
