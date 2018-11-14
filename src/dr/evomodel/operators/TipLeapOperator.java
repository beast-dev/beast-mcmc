package dr.evomodel.operators;


import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.operators.TipLeapOperatorParser;
import dr.inference.operators.*;
import dr.math.MathUtils;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Implements the Tip Leap move.
 *
 * This move picks a taxon at random from a list of taxa, find the corresponding node, and then moves the parent to any location
 * that is a certain patristic distance from its starting point (the distance is drawn from a Gaussian).
 *
 * It is always possible for the node to move up (potentially becoming the root) but the destination can't
 * be younger than the original node. All possible destinations are collected and then picked amongst
 * uniformly.
 *
 * This operator is adapted from SubtreeLeapMove
 *
 * @author Mathieu Fourment
 * @version $Id$
 */

public class TipLeapOperator extends AbstractTreeOperator implements AdaptableMCMCOperator {
    private double size;

    private final TreeModel tree;
    private final AdaptationMode mode;
    private final double targetAcceptance;

    private Taxa taxa;

    /**
     * Constructor
     *
     * @param tree   the tree
     * @param taxa   some taxa
     * @param weight the weight
     * @param size   scaling on a unit Gaussian to draw the patristic distance from
     * @param mode   coercion mode
     */
    public TipLeapOperator(TreeModel tree, Taxa taxa, double weight, double size, double targetAcceptance, AdaptationMode mode) {
        this.tree = tree;
        setWeight(weight);
        this.size = size;
        this.targetAcceptance = targetAcceptance;
        this.mode = mode;
        this.taxa = taxa;
    }



    /**
     * Do a tip leap move.
     *
     * @return the log-transformed hastings ratio
     */
    public double doOperation() {
        double logq;

        final double delta = getDelta();

        // Pick a taxon from taxa
        Taxon taxon = taxa.getTaxon(MathUtils.nextInt(taxa.getTaxonCount()));
        NodeRef node = null;
        int external = tree.getExternalNodeCount();
        for (int i = 0; i < external; i++) {
            if (tree.getNodeTaxon(tree.getExternalNode(i)).getId().equals(taxon.getId())) {
                node = tree.getExternalNode(i);
                break;
            }
        }
        if (node == null) {
            System.err.println("Taxon not found in tree: " + taxon.getId());
            System.exit(1);
        }

        // get its parent - this is the node we will prune/graft
        final NodeRef parent = tree.getParent(node);

        // get the node's sibling
        final NodeRef sibling = getOtherChild(tree, parent, node);

        // and its grand parent
        final NodeRef grandParent = tree.getParent(parent);

        final Map<NodeRef, Double> destinations = getDestinations(node, parent, sibling, delta);
        final List<NodeRef> destinationNodes = new ArrayList<NodeRef>(destinations.keySet());

        // pick uniformly from this list
        int r = MathUtils.nextInt(destinations.size());

        double forwardProbability = 1.0 / destinations.size();

        final NodeRef j = destinationNodes.get(r);
        final double newHeight = destinations.get(j);

        final NodeRef jParent = tree.getParent(j);

        if (jParent != null && newHeight > tree.getNodeHeight(jParent)) {
            throw new IllegalArgumentException("height error");
        }

        if (newHeight < tree.getNodeHeight(j)) {
            throw new IllegalArgumentException("height error");
        }

        tree.beginTreeEdit();

        if (j == parent || jParent == parent) {
            // the subtree is not actually moving but the height will change
        } else {
            if (grandParent == null) {
                // if the parent of the original node is the root then the sibling becomes
                // the root.
                tree.removeChild(parent, sibling);
                tree.setRoot(sibling);

            } else {
                // remove the parent of node by connecting its sibling to its grandparent.
                tree.removeChild(parent, sibling);
                tree.removeChild(grandParent, parent);
                tree.addChild(grandParent, sibling);
            }

            if (jParent == null) {
                // adding the node to the root of the tree
                tree.addChild(parent, j);
                tree.setRoot(parent);
            } else {
                // remove destination edge j from its parent
                tree.removeChild(jParent, j);

                // add destination edge to the parent of node
                tree.addChild(parent, j);

                // and add the parent of i as a child of the former parent of j
                tree.addChild(jParent, parent);
            }
        }
        tree.endTreeEdit();

        tree.setNodeHeight(parent, newHeight);

        if (tree.getParent(parent) != null && newHeight > tree.getNodeHeight(tree.getParent(parent))) {
            throw new IllegalArgumentException("height error");
        }

        if (newHeight < tree.getNodeHeight(node)) {
            throw new IllegalArgumentException("height error");
        }

        if (newHeight < tree.getNodeHeight(getOtherChild(tree, parent, node))) {
            throw new IllegalArgumentException("height error");
        }

        final Map<NodeRef, Double> reverseDestinations = getDestinations(node, parent, getOtherChild(tree, parent, node), delta);
        double reverseProbability = 1.0 / reverseDestinations.size();

        // hastings ratio = reverse Prob / forward Prob
        logq = Math.log(reverseProbability) - Math.log(forwardProbability);
        return logq;
    }

    private Map<NodeRef, Double> getDestinations(NodeRef node, NodeRef parent, NodeRef sibling, double delta) {

        final Map<NodeRef, Double> destinations = new HashMap<NodeRef, Double>();

        // get the parent's height
        final double height = tree.getNodeHeight(parent);

        final double heightBelow = height - delta;

        if (heightBelow > tree.getNodeHeight(node)) {
            // the destination height below the parent is compatible with the node
            // see if there are any destinations on the sibling's branch
            final List<NodeRef> edges = new ArrayList<NodeRef>();

            getIntersectingEdges(tree, sibling, heightBelow, edges);

            // add the intersecting edges and the height
            for (NodeRef n : edges) {
                destinations.put(n, heightBelow);
            }
        }

        final double heightAbove = height + delta;

        NodeRef node1 = parent;

        // walk up to root
        boolean done = false;
        while (!done) {
            NodeRef parent1 = tree.getParent(node1);

            if (parent1 != null) {
                final double height1 = tree.getNodeHeight(parent1);
                if (height1 < heightAbove) {
                    // haven't reached the height above the original height so go down
                    // the sibling subtree
                    NodeRef sibling1 = getOtherChild(tree, parent1, node1);

                    double heightBelow1 = height1 - (heightAbove - height1);

                    if (heightBelow1 > tree.getNodeHeight(node)) {

                        final List<NodeRef> edges = new ArrayList<NodeRef>();

                        getIntersectingEdges(tree, sibling1, heightBelow1, edges);

                        // add the intersecting edges and the height
                        for (NodeRef n : edges) {
                            destinations.put(n, heightBelow1);
                        }
                    }
                } else {
                    // add the current node as a destination
                    destinations.put(node1, heightAbove);
                    done = true;
                }

                node1 = parent1;
            } else {
                // node1 is the root - add it as a destination and stop loop
                destinations.put(node1, heightAbove);
                done = true;
            }
        }

        return destinations;
    }

    private double getDelta() {
        return Math.abs(MathUtils.nextGaussian() * size);
    }

    private int getIntersectingEdges(Tree tree, NodeRef node, double height, List<NodeRef> edges) {

        final NodeRef parent = tree.getParent(node);

        if (tree.getNodeHeight(parent) < height) return 0;

        if (tree.getNodeHeight(node) < height) {
            edges.add(node);
            return 1;
        }

        int count = 0;
        for (int i = 0; i < tree.getChildCount(node); i++) {
            count += getIntersectingEdges(tree, tree.getChild(node, i), height, edges);
        }
        return count;
    }

    public double getSize() {
        return size;
    }

    public void setSize(double size) {
        this.size = size;
    }

    public double getAdaptableParameter() {
        return Math.log(getSize());
    }

    public void setAdaptableParameter(double value) {
        setSize(Math.exp(value));
    }

    public double getRawParameter() {
        return getSize();
    }

    public AdaptationMode getMode() {
        return mode;
    }

    public double getMinimumAcceptanceLevel() {
        return 0.1;
    }

    public double getMaximumAcceptanceLevel() {
        return 0.4;
    }

    public double getMinimumGoodAcceptanceLevel() {
        return 0.20;
    }

    public double getMaximumGoodAcceptanceLevel() {
        return 0.30;
    }

    public final String getPerformanceSuggestion() {
        return null;
    }

    @Override
    public double getTargetAcceptanceProbability() {
        return targetAcceptance;
    }

    public String getAdaptableParameterName() {
        return "size";
    }

    public String getOperatorName() {
        return TipLeapOperatorParser.TIP_LEAP + "(" + tree.getId() + ")";
    }


}