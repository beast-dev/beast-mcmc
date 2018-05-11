/*
 * SubtreeLeapOperator.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.evomodel.operators;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.operators.SubtreeLeapOperatorParser;
import dr.inference.operators.*;
import dr.math.MathUtils;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Implements the Subtree Leap move.
 *
 * This move picks a node at random (except for the root) and then moves the parent to any location
 * that is a certain patristic distance from its starting point (the distance is drawn from a Gaussian).
 *
 * It is always possible for the node to move up (potentially becoming the root) but the destination can't
 * be younger than the original node. All possible destinations are collected and then picked amongst
 * uniformly.
 *
 * @author Andrew Rambaut
 * @version $Id$
 */

public class SubtreeLeapOperator extends AbstractTreeOperator implements CoercableMCMCOperator {

    private static final boolean DEBUG = false;
    
    private double size = 1.0;
    private double accP = 0.234;

    private final TreeModel tree;
    private final CoercionMode mode;

    private Instance lastInstance;

    /**
     * Constructor
     *
     * @param tree   the tree
     * @param weight the weight
     * @param size   scaling on a unit Gaussian to draw the patristic distance from
     * @param mode   coercion mode
     */
    public SubtreeLeapOperator(TreeModel tree, double weight, double size, double accP, CoercionMode mode) {
        this.tree = tree;
        setWeight(weight);
        this.size = size;
        this.accP = accP;
        this.mode = mode;
    }

    public TreeModel getTreeModel() {
        return tree;
    }

    /**
     * Do a subtree leap move.
     *
     * @return the log-transformed hastings ratio
     */
    public double doOperation() {

        lastInstance = drawOperation();

        return applyInstance(lastInstance);
    }

    protected Instance drawOperation() {

        final double delta = getDelta();

        final NodeRef root = tree.getRoot();

        NodeRef node;

        do {
            // choose a random node avoiding root
            node = tree.getNode(MathUtils.nextInt(tree.getNodeCount()));

        } while (node == root);

        // get its parent - this is the node we will prune/graft
        final NodeRef parent = tree.getParent(node);

        // get the node's sibling
        final NodeRef sibling = getOtherChild(tree, parent, node);

        // and its grand parent
        final NodeRef grandParent = tree.getParent(parent);

        final Map<NodeRef, Double> destinations = getDestinations(node, parent, sibling, delta);
        final List<NodeRef> destinationNodes = new ArrayList<NodeRef>(destinations.keySet());

        double forwardProbability = 1.0 / destinations.size();

        // pick uniformly from this list
        int r = MathUtils.nextInt(destinations.size());

        final NodeRef destination = destinationNodes.get(r);
        final double destinationHeight = destinations.get(destination);

        final NodeRef destinationParent = tree.getParent(destination);

        return new Instance(delta, node, parent, sibling, grandParent, destination, destinationParent, destinationHeight, forwardProbability);
    }

    protected Instance getLastInstance() {
        return lastInstance;
    }

    protected double applyInstance(Instance instance) {
        tree.beginTreeEdit();

        NodeRef parent = instance.parent < 0 ? null : tree.getNode(instance.parent);
        NodeRef grandParent = instance.grandParent < 0 ? null : tree.getNode(instance.grandParent);
        NodeRef sibling = instance.sibling < 0 ? null : tree.getNode(instance.sibling);
        NodeRef destination = instance.destination < 0 ? null : tree.getNode(instance.destination);
        NodeRef destinationParent = instance.destinationParent < 0 ? null : tree.getNode(instance.destinationParent);

        if (destination == parent || destinationParent == parent) {
            // the subtree is not actually moving but the height will change

            assert true;
        } else {
            if (grandParent == null) {
                // if the parent of the original node is the root then the sibling becomes
                // the root.
                if (DEBUG) System.err.println("STL: removing child " + sibling.getNumber() + " of root node " + parent.getNumber());

                tree.removeChild(parent, sibling);
                if (DEBUG) System.err.println("STL: setting node " + sibling.getNumber() + " as root node ");

                tree.setRoot(sibling);

            } else {
                // remove the parent of node by connecting its sibling to its grandparent.
                if (DEBUG) System.err.println("STL: removing child " + sibling.getNumber() + " of node " + parent.getNumber() );
                tree.removeChild(parent, sibling);
                if (DEBUG) System.err.println("STL: removing child " + parent.getNumber() + " of node " + grandParent.getNumber() );
                tree.removeChild(grandParent, parent);
                if (DEBUG) System.err.println("STL: setting node " + sibling.getNumber() + " as child of node " + grandParent.getNumber());
                tree.addChild(grandParent, sibling);
            }

            if (destinationParent == null) {
                // adding the node to the root of the tree
                if (DEBUG) System.err.println("STL: adding node " + destination.getNumber() + " as child of node " + parent.getNumber());
                tree.addChild(parent, destination);
                if (DEBUG) System.err.println("STL: setting node " + parent.getNumber() + " as root node ");
                tree.setRoot(parent);
            } else {
                if (DEBUG) System.err.println("STL: removing child " + destination.getNumber() + " of node " + destinationParent.getNumber() );
                // remove destination edge from its parent
                tree.removeChild(destinationParent, destination);

                if (DEBUG) System.err.println("STL: adding node " + destination.getNumber() + " as child of node " + parent.getNumber());
                // add destination edge to the parent of node
                tree.addChild(parent, destination);

                if (DEBUG) System.err.println("STL: adding node " + parent.getNumber() + " as child of node " + destinationParent.getNumber());
                // and add the parent of target node as a child of the former parent of destination
                tree.addChild(destinationParent, parent);
            }
        }
        tree.endTreeEdit();

        if (DEBUG) System.err.println("STL: setting height of node " + parent.getNumber() + " to " + instance.destinationHeight);
        tree.setNodeHeight(parent, instance.destinationHeight);

        NodeRef node = tree.getNode(instance.node);

        final Map<NodeRef, Double> reverseDestinations = getDestinations(
                node, parent, getOtherChild(tree, parent, node), instance.delta);
        double reverseProbability = 1.0 / reverseDestinations.size();

        // hastings ratio = reverse Prob / forward Prob
        double logHastingsRatio = Math.log(reverseProbability) - Math.log(instance.forwardProbability);

        return logHastingsRatio;
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

    public double getCoercableParameter() {
        return Math.log(getSize());
    }

    public void setCoercableParameter(double value) {
        setSize(Math.exp(value));
    }

    public double getRawParameter() {
        return getSize();
    }

    public CoercionMode getMode() {
        return mode;
    }

    public double getTargetAcceptanceProbability() {
        return accP;
    }


    public String getPerformanceSuggestion() {
        double prob = MCMCOperator.Utils.getAcceptanceProbability(this);
        double targetProb = getTargetAcceptanceProbability();

        double ws = OperatorUtils.optimizeWindowSize(getSize(), Double.MAX_VALUE, prob, targetProb);

        if (prob < getMinimumGoodAcceptanceLevel()) {
            return "Try decreasing size to about " + ws;
        } else if (prob > getMaximumGoodAcceptanceLevel()) {
            return "Try increasing size to about " + ws;
        } else return "";
    }

    public String getOperatorName() {
        return SubtreeLeapOperatorParser.SUBTREE_LEAP + "(" + tree.getId() + ")";
    }

    class Instance {
        final double delta;
        final int node;
        final int parent;
        final int sibling;
        final int grandParent;
        final int destination;
        final int destinationParent;
        final double destinationHeight;
        final double forwardProbability;

        public Instance(double delta, NodeRef node, NodeRef parent, NodeRef sibling, NodeRef grandParent, NodeRef destination, NodeRef destinationParent, double destinationHeight, double forwardProbability) {
            this.delta = delta;
            this.node = node == null ? -1 : node.getNumber();
            this.parent = parent == null ? -1 : parent.getNumber();
            this.sibling = sibling == null ? -1 : sibling.getNumber();
            this.grandParent = grandParent == null ? -1 : grandParent.getNumber();
            this.destination = destination == null ? -1 : destination.getNumber();
            this.destinationParent = destinationParent == null ? -1 : destinationParent.getNumber();
            this.destinationHeight = destinationHeight;
            this.forwardProbability = forwardProbability;
        }
    }
}