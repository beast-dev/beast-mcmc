package dr.evomodel.treelikelihood.thorneytreelikelihood;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;

public interface BranchLengthProvider  {
    double getBranchLength(final Tree tree, final NodeRef node);
}
