package dr.evomodel.bigfasttree;

import dr.evolution.tree.NodeRef;
import dr.evomodel.operators.AbstractTreeOperator;
import dr.evomodel.tree.TreeModel;
import dr.math.MathUtils;

import java.util.ArrayList;
import java.util.List;

/**
 * An operator that prunes a node at random and then regrafts it uniformily
 * into the same clade
 *
 * @author JT McCrone
 */
public class CladeAwareSubtreePruneRegraft extends AbstractTreeOperator {

    private final TreeModel tree;
    private final CladeNodeModel cladeModel;
    private final static String CLADE_AWARE_SUBTREE_PRUNE_REGRAFT ="CladeAwareSubtreePruneRegraft";
    private double totalDistanceTraversed = 0;
    private final static boolean DEBUG =false;
    private final int SPRperCall;


    public CladeAwareSubtreePruneRegraft(CladeNodeModel cladeModel, double weight,int SPRperCall ){
        setWeight(weight);
        this.cladeModel=cladeModel;
        this.tree = cladeModel.getTreeModel();

        this.SPRperCall = SPRperCall;

    }
    public CladeAwareSubtreePruneRegraft(CladeNodeModel cladeModel, double weight){
        this(cladeModel, weight, 1);

    }

    @Override
    public String getOperatorName() {
        return CLADE_AWARE_SUBTREE_PRUNE_REGRAFT;
    }

    @Override
    public double doOperation() {

        if(DEBUG){
            System.out.println("Before move " );
//            System.out.println(tree.toString());
        }

        CladeRef selectedClade = null; // the clade node draws can come from
        double selectedCladeSize = Double.POSITIVE_INFINITY;

        // choose a random nodes avoiding root and root's children
        // get its parent - this is the node we will prune/graft

        int chosenNodes = 0;
        while(chosenNodes<Math.min(SPRperCall,selectedCladeSize)) {
            final NodeRef root = tree.getRoot();
            NodeRef node;
            if(selectedClade==null){
                do {
                    node = tree.getNode(MathUtils.nextInt(tree.getNodeCount()));
                } while (node == root || tree.getParent(node) == root);
                //Set clade and update number of draws if needed
                selectedClade = cladeModel.getClade(node);
                selectedCladeSize = cladeModel.getNodeCount(selectedClade);
                if (cladeModel.getParent(selectedClade) == null) {
                    selectedCladeSize -= 3; // can't have root or children and this is the root clade
                }
            }else{
                //Choose the rest from the same clade so we limit updates on tree
                do {
                    node = cladeModel.getNode(selectedClade, MathUtils.nextInt((int) selectedCladeSize));
                } while (node == root || tree.getParent(node) == root);
            }

            chosenNodes++;

            if(DEBUG){
                System.out.println("setting up family for move " + chosenNodes);
            }

            final NodeRef parent = tree.getParent(node);
            final CladeRef eligibleClade = cladeModel.getClade(parent); // the clade that defined possible positions
            final NodeRef sibling = getOtherChild(tree, parent, node);
            final NodeRef grandParent = tree.getParent(parent);

            if(DEBUG){
                System.out.println("removing parent");
            }

            // prune parent out of tree
            tree.beginTreeEdit();
            // if the parent of the original node is the root then the sibling becomes
            // the root.

            // remove the parent of node by connecting its sibling to its grandparent.
            tree.removeChild(parent, sibling);
            tree.removeChild(grandParent, parent);
            tree.addChild(grandParent, sibling);
            // if parent is the root of the clade make the sibling the root of the clade for now
            if (parent == cladeModel.getRootNode(eligibleClade)) {
                cladeModel.setRootNode(eligibleClade, sibling);
            }
            totalDistanceTraversed = 0;

            if(DEBUG){
                System.out.println("getting distances");
            }

            List<Double> nodeDistances = new ArrayList<>();
            List<Double> effectiveNodeHeights = new ArrayList<>();
            List<NodeRef> nodes = new ArrayList<>();

            getDistances(tree.getNodeHeight(node), cladeModel.getRootNode(eligibleClade), eligibleClade, nodeDistances, effectiveNodeHeights, nodes);

            if(DEBUG){
                System.out.println("choosing node index");
            }

            // Pick uniform between 0 and max of traversal
            double u = MathUtils.nextDouble() * totalDistanceTraversed;

            int i = 0;
            double distanceTraveled = nodeDistances.get(i);
            while (distanceTraveled < u) {
                i++;
                distanceTraveled = nodeDistances.get(i);
            }

            if(DEBUG){
                System.out.println("choosing node");
            }

            NodeRef j = nodes.get(i);

            if(DEBUG){
                System.out.println("getting height");
            }
            double newHeight = effectiveNodeHeights.get(i) + (nodeDistances.get(i) - u);


            if (DEBUG) {
                System.out.println("node = " + node);
                System.out.println("parent= " +parent);
                System.out.println("sibling = " +sibling);
                System.out.println("newHeight: " + newHeight);
                System.out.println("destination"+ j);
            }

            NodeRef jParent = tree.getParent(j);

            if (jParent == null) {
                // parent will be the root
                tree.addChild(parent, j);
                tree.setRoot(parent);
                cladeModel.setRootNode(eligibleClade, parent);

                assert cladeModel.getParent(eligibleClade) == null;
            } else {
                // remove destination edge j from its parent
                tree.removeChild(jParent, j);

                // add destination edge to the parent of node
                tree.addChild(parent, j);

                // and add the parent of i as a child of the former parent of j
                tree.addChild(jParent, parent);

                if (cladeModel.getClade(jParent) != eligibleClade) {
                    assert cladeModel.getClade(jParent) == cladeModel.getParent(eligibleClade);
                    assert cladeModel.getRootNode(eligibleClade) == j;
                    cladeModel.setRootNode(eligibleClade, parent);
                }
            }

            tree.setNodeHeight(parent, newHeight);
            tree.endTreeEdit();
            if(DEBUG){
                System.out.println("after edit");
            }
            if (tree.getParent(parent) != null && newHeight > tree.getNodeHeight(tree.getParent(parent))) {
                throw new IllegalArgumentException("height error");
            }

            if (newHeight < tree.getNodeHeight(node)) {
                throw new IllegalArgumentException("height error: new height above node");
            }

            if (newHeight < tree.getNodeHeight(getOtherChild(tree, parent, node))) {
                assert getOtherChild(tree, parent, node) == j;
                throw new IllegalArgumentException("height error: new height above destination");
            }

            if (DEBUG) {
                System.out.println("after Move");
                System.out.println(tree.toString());
            }

        }

        return 0;
    }

    private void getDistances(double height,NodeRef node,CladeRef clade, List<Double> nodeDistances,List<Double> effectiveNodeHeights, List<NodeRef> nodes) {
        double effectiveHeight = Math.max(tree.getNodeHeight(node), height);

        if (tree.getParent(node) == null) {
            totalDistanceTraversed += 0;
        }else{
            double effectiveBranchLength = tree.getNodeHeight(tree.getParent(node)) - effectiveHeight;
            totalDistanceTraversed += effectiveBranchLength;
        }
        nodeDistances.add(totalDistanceTraversed);
        nodes.add(node);
        effectiveNodeHeights.add(effectiveHeight);

        if(cladeModel.getClade(node)==clade && tree.getNodeHeight(node)>height){
            for (int i = 0; i < tree.getChildCount(node); i++) {
                NodeRef child = tree.getChild(node, i);
                getDistances(height,child,clade,nodeDistances,effectiveNodeHeights,nodes);
            }
        }
    }

}