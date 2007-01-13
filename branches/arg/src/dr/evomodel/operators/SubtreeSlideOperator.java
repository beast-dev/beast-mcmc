/*
 * SubtreeSlideOperator.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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

import dr.evolution.tree.MutableTree;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.ARGModel;
import dr.inference.operators.*;
import dr.math.MathUtils;
import dr.xml.AttributeRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

import java.util.ArrayList;

/**
 * Implements the subtree slide move.
 *
 * @author Alexei Drummond
 * @version $Id: SubtreeSlideOperator.java,v 1.15 2005/06/14 10:40:34 rambaut Exp $
 */
public class SubtreeSlideOperator extends SimpleMCMCOperator implements CoercableMCMCOperator {

    public static final String SUBTREE_SLIDE = "subtreeSlide";
    public static final String SWAP_RATES = "swapRates";
    public static final String SWAP_TRAITS = "swapTraits";
    private ARGModel tree = null;
    private double size = 1.0;
    private boolean gaussian = false;
    private boolean swapRates;
    private boolean swapTraits;
    private int mode = CoercableMCMCOperator.DEFAULT;

    public SubtreeSlideOperator(ARGModel tree, int weight, double size, boolean gaussian, boolean swapRates, boolean swapTraits, int mode) {
        this.tree = tree;
        setWeight(weight);

        this.size = size;
        this.gaussian = gaussian;
        this.swapRates = swapRates;
        this.swapTraits = swapTraits;

        this.mode = mode;
    }

    /**
     * Do a probablistic subtree slide move.
     *
     * @return the log-transformed hastings ratio
     */
    public double doOperation() throws OperatorFailedException {

        double logq;

        NodeRef i, newParent, newChild;

        // 1. choose a random node avoiding root
        do {
            i = tree.getNode(MathUtils.nextInt(tree.getNodeCount()));
        } while (tree.getRoot() == i);
        NodeRef iP = tree.getParent(i);
        NodeRef CiP = getOtherChild(tree, iP, i);
        NodeRef PiP = tree.getParent(iP);

        // 2. choose a delta to move
        double delta = getDelta();
        double oldHeight = tree.getNodeHeight(iP);
        double newHeight = oldHeight + delta;

        // 3. if the move is up
        if (delta > 0) {

            // 3.1 if the topology will change
            if (PiP != null && tree.getNodeHeight(PiP) < newHeight) {

                // find new parent
                newParent = PiP;
                newChild = iP;
                while (tree.getNodeHeight(newParent) < newHeight) {
                    newChild = newParent;
                    newParent = tree.getParent(newParent);
                    if (newParent == null) break;
                }

                tree.beginTreeEdit();

                // 3.1.1 if creating a new root
                if (tree.isRoot(newChild)) {
                    tree.removeChild(iP, CiP);
                    tree.removeChild(PiP, iP);
                    tree.addChild(iP, newChild);
                    tree.addChild(PiP, CiP);
                    tree.setRoot(iP);
                    //System.err.println("Creating new root!");
                }
                // 3.1.2 no new root
                else {
                    tree.removeChild(iP, CiP);
                    tree.removeChild(PiP, iP);
                    tree.removeChild(newParent, newChild);
                    tree.addChild(iP, newChild);
                    tree.addChild(PiP, CiP);
                    tree.addChild(newParent, iP);
                    //System.err.println("No new root!");
                }

                tree.setNodeHeight(iP, newHeight);

                try {
                    tree.endTreeEdit();
                } catch (MutableTree.InvalidTreeException ite) {
                    throw new RuntimeException(ite.toString());
                }

                // 3.1.3 count the hypothetical sources of this destination.
                int possibleSources = intersectingEdges(tree, newChild, oldHeight, null);
                //System.out.println("possible sources = " + possibleSources);

                logq = Math.log(1.0 / (double) possibleSources);

            } else {
                // just change the node height
                tree.setNodeHeight(iP, newHeight);
                logq = 0.0;
            }
        }
        // 4 if we are sliding the subtree down.
        else {

            // 4.0 is it a valid move?
            if (tree.getNodeHeight(i) > newHeight) {
                return Double.NEGATIVE_INFINITY;
            }

            // 4.1 will the move change the topology
            if (tree.getNodeHeight(CiP) > newHeight) {

                ArrayList newChildren = new ArrayList();
                int possibleDestinations = intersectingEdges(tree, CiP, newHeight, newChildren);

                // if no valid destinations then return a failure
                if (newChildren.size() == 0) {
                    return Double.NEGATIVE_INFINITY;
                }

                // pick a random parent/child destination edge uniformly from options
                int childIndex = MathUtils.nextInt(newChildren.size());
                newChild = (NodeRef) newChildren.get(childIndex);
                newParent = tree.getParent(newChild);

                tree.beginTreeEdit();

                // 4.1.1 if iP was root
                if (tree.isRoot(iP)) {
                    // new root is CiP
                    tree.removeChild(iP, CiP);
                    tree.removeChild(newParent, newChild);
                    tree.addChild(iP, newChild);
                    tree.addChild(newParent, iP);
                    tree.setRoot(CiP);
                    //System.err.println("DOWN: Creating new root!");
                } else {
                    tree.removeChild(iP, CiP);
                    tree.removeChild(PiP, iP);
                    tree.removeChild(newParent, newChild);
                    tree.addChild(iP, newChild);
                    tree.addChild(PiP, CiP);
                    tree.addChild(newParent, iP);
                    //System.err.println("DOWN: no new root!");
                }

                tree.setNodeHeight(iP, newHeight);

                try {
                    tree.endTreeEdit();
                } catch (MutableTree.InvalidTreeException ite) {
                    throw new RuntimeException(ite.toString());
                }

                logq = Math.log((double) possibleDestinations);
            } else {
                tree.setNodeHeight(iP, newHeight);
                logq = 0.0;
            }
        }

        if (swapRates) {
            NodeRef j = tree.getNode(MathUtils.nextInt(tree.getNodeCount()));
            if (j != i) {
                double tmp = tree.getNodeRate(i);
                tree.setNodeRate(i, tree.getNodeRate(j));
                tree.setNodeRate(j, tmp);
            }

        }

        if (swapTraits) {
            NodeRef j = tree.getNode(MathUtils.nextInt(tree.getNodeCount()));
            if (j != i) {
                double tmp = tree.getNodeTrait(i);
                tree.setNodeTrait(i, tree.getNodeTrait(j));
                tree.setNodeTrait(j, tmp);
            }

        }

        if (logq == Double.NEGATIVE_INFINITY) throw new OperatorFailedException("invalid slide");
        return logq;
    }

    private double getDelta() {
        if (!gaussian) {
            return (MathUtils.nextDouble() * size) - (size / 2.0);
        } else {
            return MathUtils.nextGaussian() * size;
        }
    }

    private int intersectingEdges(Tree tree, NodeRef node, double height, ArrayList directChildren) {

        NodeRef parent = tree.getParent(node);

        if (tree.getNodeHeight(parent) < height) return 0;

        if (tree.getNodeHeight(node) < height) {
            if (directChildren != null) directChildren.add(node);
            return 1;
        }

        int count = 0;
        for (int i = 0; i < tree.getChildCount(node); i++) {
            count += intersectingEdges(tree, tree.getChild(node, i), height, directChildren);
        }
        return count;
    }

    /**
     * @return the other child of the given parent.
     */
    private NodeRef getOtherChild(Tree tree, NodeRef parent, NodeRef child) {

        if (tree.getChild(parent, 0) == child) {
            return tree.getChild(parent, 1);
        } else {
            return tree.getChild(parent, 0);
        }
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

    public int getMode() {
        return mode;
    }

    public double getTargetAcceptanceProbability() {
        return 0.234;
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
        return SUBTREE_SLIDE;
    }

    public static dr.xml.XMLObjectParser PARSER = new dr.xml.AbstractXMLObjectParser() {

        public String getParserName() {
            return SUBTREE_SLIDE;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            boolean swapRates = false;
            boolean swapTraits = false;

            int mode = CoercableMCMCOperator.DEFAULT;
            if (xo.hasAttribute(AUTO_OPTIMIZE)) {
                if (xo.getBooleanAttribute(AUTO_OPTIMIZE)) {
                    mode = CoercableMCMCOperator.COERCION_ON;
                } else {
                    mode = CoercableMCMCOperator.COERCION_OFF;
                }
            }

            if (xo.hasAttribute(SWAP_RATES)) {
                swapRates = xo.getBooleanAttribute(SWAP_RATES);
            }
            if (xo.hasAttribute(SWAP_TRAITS)) {
                swapTraits = xo.getBooleanAttribute(SWAP_TRAITS);
            }

            ARGModel treeModel = (ARGModel) xo.getChild(ARGModel.class);
            int weight = xo.getIntegerAttribute("weight");
            double size = xo.getDoubleAttribute("size");
            boolean gaussian = xo.getBooleanAttribute("gaussian");
            return new SubtreeSlideOperator(treeModel, weight, size, gaussian, swapRates, swapTraits, mode);
        }

        public String getParserDescription() {
            return "An operator that slides a subtree.";
        }

        public Class getReturnType() {
            return SubtreeSlideOperator.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                AttributeRule.newIntegerRule("weight"),
                AttributeRule.newDoubleRule("size"),
                AttributeRule.newBooleanRule("gaussian"),
                AttributeRule.newBooleanRule(SWAP_RATES, true),
                AttributeRule.newBooleanRule(SWAP_TRAITS, true),
                AttributeRule.newBooleanRule(AUTO_OPTIMIZE, true),
//			new ElementRule(ARGModel.class)
        };
	};

}
