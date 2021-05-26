package dr.evomodel.treelikelihood.thorneytreelikelihood;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.operators.AbstractTreeOperator;
import dr.evomodel.tree.TreeModel;
import dr.math.MathUtils;

import java.util.ArrayList;
import java.util.List;

/**
 * An operator that prunes a node at random and then regrafts it a position choosen uniformily from valid branches and heights
 *
 * @author JT McCrone
 */

public class UniformSubtreePruneRegraft extends AbstractTreeOperator implements ConstrainableTreeOperator{

    private final TreeModel tree;
    private final static String UNIFORM_SUBTREE_PRUNE_REGRAFT = "UniformSubtreePruneRegraft";
    private double totalDistanceTraversed = 0;
    private final static boolean DEBUG = false;

    public UniformSubtreePruneRegraft(TreeModel tree, double weight) {
        setWeight(weight);
       this.tree = tree;

    }

    @Override
    public String getOperatorName() {
        return UNIFORM_SUBTREE_PRUNE_REGRAFT;
    }

    public double doOperation() {
        return doOperation(tree);
    }

    @Override
    public double doOperation(TreeModel tree) {


        if (DEBUG) {
            System.out.println("Before move ");
            System.out.println(tree.toString());
        }

        // choose a random nodes avoiding root and root's children
        // get its parent - this is the node we will prune/graft

        final NodeRef root = tree.getRoot();
        NodeRef node;

        do {
            node = tree.getNode(MathUtils.nextInt(tree.getNodeCount()));


        } while (node == root || tree.getParent(node) == root);

        final NodeRef parent = tree.getParent(node);

        final NodeRef sibling = getOtherChild(tree, parent, node);
        final NodeRef grandParent = tree.getParent(parent);


        if (DEBUG) {
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
        // should now be updated by clade model
//            if (parent == cladeModel.getRootNode(eligibleClade)) {
//                cladeModel.setRootNode(eligibleClade, sibling);
//            }
        totalDistanceTraversed = 0;

        if (DEBUG) {
            System.out.println("getting distances");
        }

        List<Double> nodeDistances = new ArrayList<>();
        List<Double> effectiveNodeHeights = new ArrayList<>();
        List<NodeRef> nodes = new ArrayList<>();

        getDistances(tree, tree.getNodeHeight(node), tree.getRoot(), nodeDistances, effectiveNodeHeights, nodes);

        if (DEBUG) {
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

        if (DEBUG) {
            System.out.println("choosing node");
        }

        NodeRef j = nodes.get(i);

        if (DEBUG) {
            System.out.println("getting height");
        }
        double newHeight = effectiveNodeHeights.get(i) + (nodeDistances.get(i) - u);


        if (DEBUG) {
            System.out.println("node = " + node.getNumber());
            System.out.println("parent= " + parent.getNumber());
            System.out.println("sibling = " + sibling.getNumber());
            System.out.println("newHeight = " + newHeight);
            System.out.println("destination = " + j.getNumber());
        }

        NodeRef jParent = tree.getParent(j);

        if (jParent == null) {
            // parent will be the root
            tree.addChild(parent, j);
            tree.setRoot(parent);
//            cladeModel.setRootNode(eligibleClade, parent);

        } else {
            // remove destination edge j from its parent
            tree.removeChild(jParent, j);

            // add destination edge to the parent of node
            tree.addChild(parent, j);

            // and add the parent of i as a child of the former parent of j
            tree.addChild(jParent, parent);

//            if (cladeModel.getClade(jParent) != eligibleClade) {
//                assert cladeModel.getClade(jParent) == cladeModel.getParent(eligibleClade);
//                assert cladeModel.getRootNode(eligibleClade) == j;
//                cladeModel.setRootNode(eligibleClade, parent);
//            }
        }

        tree.setNodeHeight(parent, newHeight);
        tree.endTreeEdit();
        if (DEBUG) {
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

        return 0;
    }

    private void getDistances(Tree tree, double height, NodeRef node, List<Double> nodeDistances, List<Double> effectiveNodeHeights, List<NodeRef> nodes) {
        double effectiveHeight = Math.max(tree.getNodeHeight(node), height);

        if (tree.getParent(node) == null) {
            totalDistanceTraversed += 0;
        } else {
            double effectiveBranchLength = tree.getNodeHeight(tree.getParent(node)) - effectiveHeight;
            totalDistanceTraversed += effectiveBranchLength;
        }
        nodeDistances.add(totalDistanceTraversed);
        nodes.add(node);
        effectiveNodeHeights.add(effectiveHeight);

        boolean goOn =  tree.getNodeHeight(node) > height;

        if (goOn) {
            for (int i = 0; i < tree.getChildCount(node); i++) {
                NodeRef child = tree.getChild(node, i);
                getDistances(tree, height, child,  nodeDistances, effectiveNodeHeights, nodes);
            }
        }
    }

}