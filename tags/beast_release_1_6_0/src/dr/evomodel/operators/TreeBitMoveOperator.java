package dr.evomodel.operators;

import dr.evolution.tree.NodeRef;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;

import java.util.ArrayList;
import java.util.List;

/**
 * A generic operator swaps a randomly selected rate change from parent to offspring or vice versa.
 *
 * @author Alexei Drummond
 * @version $Id$
 */
public class TreeBitMoveOperator extends SimpleMCMCOperator {

    public TreeBitMoveOperator(TreeModel tree, String t1, String t2, double weight) {
        this.tree = tree;
        this.indicatorTrait = t1;
        this.trait2 = t2;

        if (indicatorTrait == null) indicatorTrait = "trait";

        setWeight(weight);
    }

    /**
     * Pick a parent-child node pair involving a single rate change and swap the rate change location
     * and corresponding rate parameters.
     */
    public final double doOperation() throws OperatorFailedException {

        NodeRef root = tree.getRoot();

        // 1. collect nodes that form a pair with parent such that
        // one of them has a one and one has a zero
        List<NodeRef> candidates = new ArrayList<NodeRef>();

        for (int i = 0; i < tree.getNodeCount(); i++) {

            NodeRef node = tree.getNode(i);
            if (node != root && tree.getParent(node) != root) {

                NodeRef parent = tree.getParent(node);

                int sum = rateChange(tree, node) + rateChange(tree, parent);

                if (sum == 1) candidates.add(node);
            }
        }

        if (candidates.size() == 0) throw new OperatorFailedException("No suitable pairs!");

        NodeRef node = candidates.get(MathUtils.nextInt(candidates.size()));
        NodeRef parent = tree.getParent(node);

        double nodeTrait, parentTrait;
        double nodeRate, parentRate;

        nodeTrait = tree.getNodeTrait(node, indicatorTrait);
        parentTrait = tree.getNodeTrait(parent, indicatorTrait);

        tree.setNodeTrait(node, indicatorTrait, parentTrait);
        tree.setNodeTrait(parent, indicatorTrait, nodeTrait);

        if (trait2 != null) {
            nodeTrait = tree.getNodeTrait(node, trait2);
            parentTrait = tree.getNodeTrait(parent, trait2);

            tree.setNodeTrait(node, trait2, parentTrait);
            tree.setNodeTrait(parent, trait2, nodeTrait);
        } else {
            nodeRate = tree.getNodeRate(node);
            parentRate = tree.getNodeRate(parent);

            tree.setNodeRate(node, parentRate);
            tree.setNodeRate(parent, nodeRate);
        }

        return 0.0;
    }

    public final int rateChange(TreeModel tree, NodeRef node) {
        return (int) Math.round(tree.getNodeTrait(node, indicatorTrait));
    }

    // Interface MCMCOperator
    public final String getOperatorName() {
        return "treeBitMove()";
    }

    public final String getPerformanceSuggestion() {

        return "no performance suggestion";
    }

    public String toString() {
        return getOperatorName();
    }
    // Private instance variables

    private TreeModel tree;
    private String indicatorTrait;
    private String trait2;
}
