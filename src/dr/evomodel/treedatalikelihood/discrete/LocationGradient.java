package dr.evomodel.treedatalikelihood.discrete;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.inference.model.Parameter;

/**
 * @author Marc A. Suchard
 */
public class LocationGradient extends HyperParameterBranchRateGradient {

    public LocationGradient(String traitName, TreeDataLikelihood treeDataLikelihood,
                            BeagleDataLikelihoodDelegate likelihoodDelegate,
                            Parameter locationScaleParameter, boolean useHessian) {

        super(traitName, treeDataLikelihood, likelihoodDelegate, locationScaleParameter, useHessian);
    }

    @Override
    double getDifferential(Tree tree, NodeRef node) {
        double rate = branchRateModel.getBranchRate(tree, node);
        return locationScaleTransform.expLocationDifferential(rate); // TODO Move function below into here?
    }
}
