package dr.evomodel.branchratemodel;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;

public interface DifferentiableBranchRatesFullMethods extends DifferentiableBranchRates {
    Tree getTree();

    double getUntransformedBranchRate(Tree tree, NodeRef node);

}
