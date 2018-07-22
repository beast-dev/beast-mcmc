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

        double[] design = fixedEffects.getDesignVector(tree, node); // TODO Use

        return new double[] {
                locationScaleTransform.expLocationDifferential(rate)  // TODO Move function below into here?
        };
    }
}
