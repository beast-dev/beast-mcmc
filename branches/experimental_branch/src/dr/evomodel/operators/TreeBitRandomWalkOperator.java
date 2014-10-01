package dr.evomodel.operators;

import dr.evolution.tree.NodeRef;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.operators.TreeBitRandomWalkOperatorParser;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;

import java.util.ArrayList;
import java.util.List;

/**
 * A generic operator that randomly selects a 1-bit in the tree and choose a new
 * location nearby, by doing k random steps on the tree with equal weights on the
 * parent and two children for each random step. The 1 bit is swapped with the bit
 * at the new location, optionally the associated variable values are swapped as well.
 *
 * @author Alexei Drummond
 * @version $Id$
 */
public class TreeBitRandomWalkOperator extends SimpleMCMCOperator {

    public TreeBitRandomWalkOperator(TreeModel tree, String t1, String t2, double weight, int k, boolean swapTrait2) {
        this.tree = tree;
        this.indicatorTrait = t1;
        this.trait2 = t2;

        if (indicatorTrait == null) indicatorTrait = "trait";
        this.k = k;

        this.swapTrait2 = swapTrait2;

        setWeight(weight);
    }

    /**
     * Pick a parent-child node pair involving a single rate change and swap the rate change location
     * and corresponding rate parameters.
     */
    public final double doOperation() throws OperatorFailedException {

        // 1. collect nodes that form a pair with parent such that
        // one of them has a one and one has a zero
        List<NodeRef> candidates = new ArrayList<NodeRef>();

        for (int i = 0; i < tree.getNodeCount(); i++) {
            NodeRef node = tree.getNode(i);
            if (tree.getNodeTrait(node, indicatorTrait) == 1.0) candidates.add(node);
        }

        if (candidates.size() == 0) throw new OperatorFailedException("No suitable bits!");

        NodeRef node = candidates.get(MathUtils.nextInt(candidates.size()));

        NodeRef newNode = node;

        for (int i = 0; i < k; i++) {
            int randomNode = MathUtils.nextInt(3);
            if (randomNode < 2) {
                if (!tree.isExternal(newNode)) {
                    newNode = tree.getChild(newNode, randomNode);
                }
            } else if (!tree.isRoot(newNode)) {
                newNode = tree.getParent(newNode);
            }
        }

        // this shortcut avoids unnecessary likelihood calculations
        if (node == newNode) throw new OperatorFailedException("Moving to same node!");

        double nodeTrait, newTrait;
        double nodeRate, newRate;

        nodeTrait = tree.getNodeTrait(node, indicatorTrait);
        newTrait = tree.getNodeTrait(newNode, indicatorTrait);

        tree.setNodeTrait(node, indicatorTrait, newTrait);
        tree.setNodeTrait(newNode, indicatorTrait, nodeTrait);

        if (swapTrait2) {
            if (trait2 != null) {
                nodeTrait = tree.getNodeTrait(node, trait2);
                newTrait = tree.getNodeTrait(newNode, trait2);

                tree.setNodeTrait(node, trait2, newTrait);
                tree.setNodeTrait(newNode, trait2, nodeTrait);
            } else {
                nodeRate = tree.getNodeRate(node);
                newRate = tree.getNodeRate(newNode);

                tree.setNodeRate(node, newRate);
                tree.setNodeRate(newNode, nodeRate);
            }
        }

        return 0.0;
    }

    public final int rateChange(TreeModel tree, NodeRef node) {
        return (int) Math.round(tree.getNodeTrait(node, indicatorTrait));
    }

    // Interface MCMCOperator
    public final String getOperatorName() {
        return TreeBitRandomWalkOperatorParser.BIT_RANDOM_WALK_OPERATOR + "(" + indicatorTrait + ")";
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
    private int k;
    private boolean swapTrait2 = true;
}
