package dr.evomodel.treedatalikelihood.discrete;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.inference.model.Parameter;

/**
 * @author Marc A. Suchard
 */
public class ScaleGradient extends HyperParameterBranchRateGradient {

    public ScaleGradient(String traitName, TreeDataLikelihood treeDataLikelihood,
                         BeagleDataLikelihoodDelegate likelihoodDelegate,
                         Parameter locationScaleParameter, boolean useHessian) {

        super(traitName, treeDataLikelihood, likelihoodDelegate, locationScaleParameter, useHessian);
    }

    @Override
    double[] getDifferential(Tree tree, NodeRef node) {
        double rate = branchRateModel.getBranchRate(tree, node);

        // TODO Can out out of this class (I think), if we provide both transform() and inverse() here.

        double tmp = (Math.log(rate / locationScaleTransform.getLocation(tree, node)) - locationScaleTransform.getTransformMu())
                /(locationScaleTransform.getTransformSigma() * locationScaleTransform.getTransformSigma()) - 1.0;

        return new double[] {tmp * rate * locationScaleTransform.getScale(tree, node) / (1.0 + locationScaleTransform.getScale(tree, node) * locationScaleTransform.getScale(tree, node))};
    }
}
