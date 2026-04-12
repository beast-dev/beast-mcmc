package dr.evomodel.treedatalikelihood.continuous;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.continuous.MultivariateElasticModel;

public interface OUDelegateContext {

    Tree getTree();

    MultivariateElasticModel getElasticModel();

    MultivariateDiffusionModel getDiffusionModel();

    int getDim();

    boolean hasDiagonalActualization();

    boolean hasBlockDiagActualization();

    boolean isSymmetric();

    double[] getBasisD();

    double[] getBasisRotations();

    double[] getOptimalValues(NodeRef node);
}
