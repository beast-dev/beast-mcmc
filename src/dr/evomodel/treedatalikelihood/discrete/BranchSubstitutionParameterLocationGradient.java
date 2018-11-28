package dr.evomodel.treedatalikelihood.discrete;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.BranchSpecificFixedEffects;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.inference.model.BranchParameter;
import dr.inference.model.Parameter;


public class BranchSubstitutionParameterLocationGradient extends HyperParameterBranchSubstitutionParameterGradient {

    private final BranchSpecificFixedEffects fixedEffects;

    public BranchSubstitutionParameterLocationGradient(String traitName,
                                                       TreeDataLikelihood treeDataLikelihood,
                                                       BeagleDataLikelihoodDelegate likelihoodDelegate,
                                                       Parameter branchSubstitutionParameter,
                                                       BranchParameter branchParameter,
                                                       boolean useHessian,
                                                       BranchSpecificFixedEffects fixedEffects) {
        super(traitName, treeDataLikelihood, likelihoodDelegate, branchSubstitutionParameter, branchParameter,
                fixedEffects.getFixedEffectsParameter(), useHessian);

        this.fixedEffects = fixedEffects;
    }

    @Override
    double[] getDifferential(Tree tree, NodeRef node) {
        double rate = branchParameter.getParameterValue(node.getNumber());
        double[] results = fixedEffects.getDifferential(rate, tree, node);

        return results;
    }

    @Override
    public String getReport() {
        return getReport(hyperParameter);
    }
}
