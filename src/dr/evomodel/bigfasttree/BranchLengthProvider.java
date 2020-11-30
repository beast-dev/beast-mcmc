package dr.evomodel.bigfasttree;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.inference.model.Model;

public interface BranchLengthProvider  {
    double getBranchLength(final Tree tree, final NodeRef node);
}
