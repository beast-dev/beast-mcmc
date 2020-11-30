package dr.evomodel.bigfasttree;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.inference.model.*;



//TODO cache branch lengths
public class ConstrainedTreeBranchLengthProvider  implements BranchLengthProvider {
    public static final String CONSTRAINED_TREE_BRANCHLENGTH_PROVIDER = "ConstrainedTreeBranchLengthProvider";
    public ConstrainedTreeBranchLengthProvider(CladeNodeModel cladeNodeModel){
//        super(CONSTRAINED_TREE_BRANCHLENGTH_PROVIDER);

        this.cladeNodeModel=cladeNodeModel;
//        addModel(cladeNodeModel);
        // set up branch lengths
        Tree dataTree = cladeNodeModel.getCladeTree();
        externalBranchLengths = new double[dataTree.getExternalNodeCount()];
        cladeBranchLengths = new double[cladeNodeModel.getCladeCount()];



        for (int i = 0; i < dataTree.getExternalNodeCount(); i++) {
            double x = dataTree.getBranchLength(dataTree.getExternalNode(i));
            externalBranchLengths[i] = x;
        }
        for (int i = 0; i < cladeNodeModel.getCladeCount(); i++) {
            CladeRef clade = cladeNodeModel.getClade(i);
            // need to get node in dataTree
            NodeRef node = dataTree.getInternalNode(clade.getNumber());
            double x =  dataTree.getBranchLength(node);
            cladeBranchLengths[i] = x;
        }
    }


    @Override
    public double getBranchLength(Tree tree, NodeRef node) {
        assert tree== cladeNodeModel.getTreeModel();
        if(tree.isExternal(node)){
            return externalBranchLengths[node.getNumber()]; //node.number==external node number since all the tips come first
        }
        CladeRef clade = cladeNodeModel.getClade(node);
        if (cladeNodeModel.getRootNode(clade) == node) {
            return cladeBranchLengths[clade.getNumber()];
        }
        return 0d;
    }

    private final double[] cladeBranchLengths;
    private final double[] externalBranchLengths;
    private final CladeNodeModel cladeNodeModel;
}
