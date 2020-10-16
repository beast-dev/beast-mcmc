package dr.evomodel.bigfasttree;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evomodel.tree.TreeModel;

import java.util.BitSet;
import java.util.HashMap;
import java.util.Map;

//TODO use a cladeNodeModel to assign a length to each clade and external node

public class ConstrainedTreeBranchLengthProvider implements BranchLengthProvider{

    public ConstrainedTreeBranchLengthProvider(CladeNodeModel cladeNodeModel){
        this.cladeNodeModel=cladeNodeModel;
        // set up branch lengths
        Tree dataTree = cladeNodeModel.getCladeTree();
        branchLengths = new double[dataTree.getNodeCount()];

        for (int i = 0; i < dataTree.getNodeCount(); i++) {
            double x = dataTree.getBranchLength(dataTree.getNode(i));
            branchLengths[i] = x;
        }
    }

    private final double[] branchLengths;
    private final CladeNodeModel cladeNodeModel;

    @Override
    public double getBranchLength(Tree tree, NodeRef node) {
        assert tree== cladeNodeModel.getTreeModel();
        if(tree.isExternal(node)){
            return branchLengths[node.getNumber()];
        }
        CladeRef clade = cladeNodeModel.getClade(node);
        if (cladeNodeModel.getRootNode(clade) == node) {
            return branchLengths[tree.getExternalNodeCount()+ clade.getNumber()];
        }
        return 0d;

    }
}
