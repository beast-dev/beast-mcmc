package dr.evomodel.treedatalikelihood.discrete;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.loggers.Loggable;
import dr.inference.model.Parameter;
import dr.xml.Reportable;

public class LocalBranchRateGradientForDiscreteTrait extends DiscreteTraitBranchRateGradient
        implements GradientWrtParameterProvider, Reportable, Loggable {

    public LocalBranchRateGradientForDiscreteTrait(String traitName,
                                                   TreeDataLikelihood treeDataLikelihood,
                                                   BeagleDataLikelihoodDelegate likelihoodDelegate,
                                                   Parameter rateParameter,
                                                   boolean useHessian) {
        super(traitName, treeDataLikelihood, likelihoodDelegate, rateParameter, useHessian);
    }

    @Override
    public double[] getGradientLogDensity() {
        double[] nodeGradient = super.getGradientLogDensity();
        double[] result = new double[tree.getNodeCount() - 1];

        for (int i = 0; i < tree.getNodeCount(); i++) {
            NodeRef node = tree.getNode(i);
            if(!tree.isRoot(node)){
                final int parameterIndex = branchRateModel.getParameterIndexFromNode(node);
                result[parameterIndex] = getSubTreeGradient(tree, node, nodeGradient);
            }
        }
        return result;
    }

    private double getSubTreeGradient(Tree tree, NodeRef node, final double[] nodeGradient) {
        final double chainGradient =  branchRateModel.getBranchRateDifferential(tree, node);
        final double branchGradient = nodeGradient[branchRateModel.getParameterIndexFromNode(node)];
        if (tree.isExternal(node)) {
            return chainGradient * branchGradient;
        } else {
            double sum = chainGradient * branchGradient;
            for (int i = 0; i < tree.getChildCount(node); i++) {
                sum += getSubTreeGradient(tree, tree.getChild(node, i), nodeGradient);
            }
            return sum;
        }
    }
}
