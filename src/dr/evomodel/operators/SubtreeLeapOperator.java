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
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.operators.SubtreeLeapOperatorParser;
import dr.evomodelxml.operators.TipLeapOperatorParser;
import dr.inference.distribution.CauchyDistribution;
import dr.inference.operators.AdaptableMCMCOperator;
import dr.inference.operators.AdaptationMode;
import dr.math.MathUtils;
import dr.math.distributions.Distribution;

import java.util.ArrayList;
import java.util.LinkedHashMap;
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
 * @author Luiz Max Carvalho
 * @author Mathieu Fourment
 * @version $Id$
 */
public class SubtreeLeapOperator extends AbstractAdaptableTreeOperator {

    public enum DistanceKernelType {
        NORMAL("normal") {
            @Override
            double getDelta(double size) {
                return Math.abs(MathUtils.nextGaussian() * size);
            }
        },
        CAUCHY("Cauchy") {
            @Override
            double getDelta(double size) {
                Distribution distK = new CauchyDistribution(0, size);
                double u = MathUtils.nextDouble();
                return Math.abs(distK.quantile(u));
            }
        };

        DistanceKernelType(String name) {
            this.name = name;
        }

        @Override
        public String toString() {
            return name;
        }

        String name;

        abstract double getDelta(double size);
    }

    private double size;

    private final TreeModel tree;
    private final DistanceKernelType distanceKernel;

    private final List<NodeRef> tips;

    /**
     * Constructor
     *
     * @param tree   the tree
     * @param weight the weight
     * @param size   scaling on a unit Gaussian to draw the patristic distance from
     * @param targetAcceptance the desired acceptance probability
     * @param distanceKernel the distribution from which to draw the patristic distance 
     * @param mode   coercion mode
     */
    public SubtreeLeapOperator(TreeModel tree, double weight, double size, DistanceKernelType distanceKernel, AdaptationMode mode, double targetAcceptance) {
        super(mode, targetAcceptance);

        this.tree = tree;
        setWeight(weight);
        this.size = size;
        this.distanceKernel = distanceKernel;
        this.tips = null;
    }

    /**
     * Constructor that takes a taxon set to pick from for the move.
     *
     * @param tree   the tree
     * @param taxa   some taxa
     * @param weight the weight
     * @param size   scaling on a unit Gaussian to draw the patristic distance from
     * @param mode   coercion mode
     */
    public SubtreeLeapOperator(TreeModel tree, TaxonList taxa, double weight, double size, DistanceKernelType distanceKernel, AdaptationMode mode, double targetAcceptance) {
        super(mode, targetAcceptance);

        this.tree = tree;
        setWeight(weight);
        this.size = size;
        this.distanceKernel = distanceKernel;
        this.tips = new ArrayList<NodeRef>();

        for (Taxon taxon : taxa) {
            boolean found = false;
            for (int i = 0; i < tree.getExternalNodeCount(); i++) {
                NodeRef tip = tree.getExternalNode(i);
                if (tree.getNodeTaxon(tip).equals(taxon)) {
                    tips.add(tip);
                    found = true;
                    break;
                }
            }

            if (!found) {
                throw new IllegalArgumentException("Taxon, " + taxon.getId() + ", not found in tree with id " + tree.getId());
            }
        }
    }


    /**
     * Do a subtree leap move.
     *
     * @return the log-transformed hastings ratio
     */
    public double doOperation() {
        double logq;

        final double delta = distanceKernel.getDelta(size);

        final NodeRef root = tree.getRoot();

        NodeRef node;

        if (tips == null) {
            // Pick a node (but not the root)
            do {
                // choose a random node avoiding root
                node = tree.getNode(MathUtils.nextInt(tree.getNodeCount()));

            } while (node == root);
        } else {
            // Pick a tip from the specified set of tips.
            node = tips.get(MathUtils.nextInt(tips.size()));
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

        final Map<NodeRef, Double> destinations = new LinkedHashMap<NodeRef, Double>();

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

    @Override
    protected void setAdaptableParameterValue(double value) {
        setSize(Math.exp(value));
    }

    @Override
    protected double getAdaptableParameterValue() {
        return Math.log(getSize());
    }

    @Override
    public double getRawParameter() {
        return getSize();
    }

    public String getAdaptableParameterName() {
        return "size";
    }

    public String getOperatorName() {
        if (tips == null) {
            return SubtreeLeapOperatorParser.SUBTREE_LEAP + "(" + tree.getId() + ")";
        } else {
            return TipLeapOperatorParser.TIP_LEAP + "(" + tree.getId() + ")";
        }
    }
}
