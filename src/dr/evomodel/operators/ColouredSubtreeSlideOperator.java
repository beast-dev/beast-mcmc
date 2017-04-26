/*
 * ColouredSubtreeSlideOperator.java
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
import dr.evomodel.coalescent.structure.ColourSamplerModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.CoercableMCMCOperator;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorUtils;
import dr.math.MathUtils;
import dr.xml.*;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

import java.util.ArrayList;
import java.util.List;

/**
 * Implements the subtree slide move.
 *
 * @author Alexei Drummond
 * @version $Id: ColouredSubtreeSlideOperator.java,v 1.4 2006/09/11 09:33:01 gerton Exp $
 */
// Cleaning out untouched stuff. Can be resurrected if needed
@Deprecated
public class ColouredSubtreeSlideOperator extends AbstractTreeOperator implements CoercableMCMCOperator {

    public static final String SUBTREE_SLIDE = "colouredSubtreeSlide";
    public static final String SWAP_RATES = "swapRates";
    public static final String SWAP_TRAITS = "swapTraits";
    private TreeModel tree = null;
    private double size = 1.0;
    private boolean gaussian = false;
    private boolean swapRates;
    private boolean swapTraits;
    private CoercionMode mode = CoercionMode.DEFAULT;

    private ColourSamplerModel colouringModel;

    public ColouredSubtreeSlideOperator(TreeModel tree, ColourSamplerModel colouringModel, double weight, double size, boolean gaussian, boolean swapRates, boolean swapTraits, CoercionMode mode) {
        this.tree = tree;
        setWeight(weight);

        this.size = size;
        this.gaussian = gaussian;
        this.swapRates = swapRates;
        this.swapTraits = swapTraits;

        this.mode = mode;

        this.colouringModel = colouringModel;

    }

    /**
     * Do a probablistic subtree slide move.
     *
     * @return the log-transformed hastings ratio
     */
    public double doOperation() {

        double logP = colouringModel.getTreeColouringWithProbability().getLogProbabilityDensity();

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
                    tree.removeChild(PiP, iP);  // remove iP from in between PiP and CiP
                    tree.addChild(iP, newChild);                           // add new edge from iP to newChild (old root)
                    tree.addChild(PiP, CiP);                               // formerly two edges
                    tree.setRoot(iP);
                    //System.err.println("Creating new root!");
                }
                // 3.1.2 no new root
                else {
                    tree.removeChild(iP, CiP);
                    tree.removeChild(PiP, iP);  // remove iP from in between PiP and CiP
                    tree.removeChild(newParent, newChild);                 // split edge: first remove old one
                    tree.addChild(iP, newChild);                           // split edge: put iP in middle
                    tree.addChild(PiP, CiP);                               // formerly two edges
                    tree.addChild(newParent, iP);                          // split edge: put iP in middle
                    //System.err.println("No new root!");
                }

                tree.setNodeHeight(iP, newHeight);

                tree.endTreeEdit();

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

                List<NodeRef> newChildren = new ArrayList<NodeRef>();
                int possibleDestinations = intersectingEdges(tree, CiP, newHeight, newChildren);

                // if no valid destinations then return a failure
                if (newChildren.size() == 0) {
                    return Double.NEGATIVE_INFINITY;
                }

                // pick a random parent/child destination edge uniformly from options
                int childIndex = MathUtils.nextInt(newChildren.size());
                newChild = newChildren.get(childIndex);
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

                tree.endTreeEdit();

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
                double tmp = tree.getNodeTrait(i, "trait");
                tree.setNodeTrait(i, "trait", tree.getNodeTrait(j, "trait"));
                tree.setNodeTrait(j, "trait", tmp);
            }

        }

        if (logq == Double.NEGATIVE_INFINITY) throw new RuntimeException("invalid slide");

        colouringModel.resample();

        double logQ = colouringModel.getTreeColouringWithProbability().getLogProbabilityDensity();

        return logq + logP - logQ;
    }

    private double getDelta() {
        if (!gaussian) {
            return (MathUtils.nextDouble() * size) - (size / 2.0);
        } else {
            return MathUtils.nextGaussian() * size;
        }
    }

    private int intersectingEdges(Tree tree, NodeRef node, double height, List<NodeRef> directChildren) {

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
        return 0.234;
    }


    public String getPerformanceSuggestion() {
        double prob = Utils.getAcceptanceProbability(this);
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

    public Element createOperatorElement(Document d) {
        return d.createElement(SUBTREE_SLIDE);
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return SUBTREE_SLIDE;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            boolean swapRates = xo.getAttribute(SWAP_RATES, false);
            boolean swapTraits = xo.getAttribute(SWAP_TRAITS, false);

            CoercionMode mode = CoercionMode.parseMode(xo);
            TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
            ColourSamplerModel colourSamplerModel = (ColourSamplerModel) xo.getChild(ColourSamplerModel.class);
            double weight = xo.getDoubleAttribute("weight");
            double size = xo.getDoubleAttribute("size");
            boolean gaussian = xo.getBooleanAttribute("gaussian");
            return new ColouredSubtreeSlideOperator(treeModel, colourSamplerModel, weight, size, gaussian, swapRates, swapTraits, mode);
        }

        public String getParserDescription() {
            return "An operator that slides a subtree.";
        }

        public Class getReturnType() {
            return ColouredSubtreeSlideOperator.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                AttributeRule.newDoubleRule("weight"),
                AttributeRule.newDoubleRule("size"),
                AttributeRule.newBooleanRule("gaussian"),
                AttributeRule.newBooleanRule(SWAP_RATES, true),
                AttributeRule.newBooleanRule(SWAP_TRAITS, true),
                AttributeRule.newBooleanRule(AUTO_OPTIMIZE, true),
                new ElementRule(TreeModel.class),
                new ElementRule(ColourSamplerModel.class)
        };
    };

}
