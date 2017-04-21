/*
 * TransmissionSubtreeSlideA.java
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

package dr.evomodel.epidemiology.casetocase.operators;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.epidemiology.casetocase.AbstractCase;
import dr.evomodel.epidemiology.casetocase.BranchMapModel;
import dr.evomodel.epidemiology.casetocase.CaseToCaseTreeLikelihood;
import dr.evomodel.operators.AbstractTreeOperator;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.*;
import dr.math.MathUtils;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 * Implements the subtree slide move where the slide does not affect the transmission tree.
 *
 * @author Matthew Hall
 *
 */
public class TransmissionSubtreeSlideA extends AbstractTreeOperator implements CoercableMCMCOperator {

    private CaseToCaseTreeLikelihood c2cLikelihood;
    private TreeModel tree;
    private double size = 1.0;
    private boolean gaussian = false;
    private final boolean swapInRandomRate;
    private final boolean swapInRandomTrait;
    private final boolean resampleInfectionTimes;
    private CoercionMode mode = CoercionMode.DEFAULT;
    private static final boolean DEBUG = false;
    public static final String TRANSMISSION_SUBTREE_SLIDE_A = "transmissionSubtreeSlideA";
    public static final String SWAP_RATES = "swapInRandomRate";
    public static final String SWAP_TRAITS = "swapInRandomTrait";
    public static final String TARGET_ACCEPTANCE = "targetAcceptance";


    public TransmissionSubtreeSlideA(CaseToCaseTreeLikelihood c2cLikelihood, double weight, double size,
                                     boolean gaussian, boolean swapRates, boolean swapTraits,  CoercionMode mode,
                                     boolean resampleInfectionTimes) {
        this.c2cLikelihood = c2cLikelihood;
        tree = c2cLikelihood.getTreeModel();
        setWeight(weight);

        if (size == 0.0) {
            double b = 0.0;
            for (int k = 0; k < tree.getNodeCount(); ++k) {
                b += tree.getBranchLength(tree.getNode(k));
            }
            size = b / (2 * tree.getNodeCount());
        }

        this.size = size;
        this.gaussian = gaussian;
        this.swapInRandomRate = swapRates;
        this.swapInRandomTrait = swapTraits;

        this.resampleInfectionTimes = resampleInfectionTimes;

        this.mode = mode;
    }

    /**
     * Do a probabilistic subtree slide move.
     *
     * @return the log-transformed hastings ratio
     */
    public double doOperation() {

        if(DEBUG){
            c2cLikelihood.outputTreeToFile("beforeTSSA.nex", false);
        }

        BranchMapModel branchMap = c2cLikelihood.getBranchMap();

        double logq = 0;

        NodeRef i;

        // 1. choose a random eligible node

        ArrayList<NodeRef> eligibleNodes = getEligibleNodes(tree, branchMap);

        i = eligibleNodes.get(MathUtils.nextInt(eligibleNodes.size()));

        double eligibleNodeCount = eligibleNodes.size();


        final NodeRef iP = tree.getParent(i);
        final NodeRef CiP = getOtherChild(tree, iP, i);
        final NodeRef PiP = tree.getParent(iP);

        // 2. choose a delta to move
        final double delta = getDelta();
        final double oldHeight = tree.getNodeHeight(iP);
        final double newHeight = oldHeight + delta;

        AbstractCase iCase = branchMap.get(i.getNumber());
        AbstractCase iPCase = branchMap.get(iP.getNumber());
        AbstractCase CiPCase = branchMap.get(CiP.getNumber());
        AbstractCase PiPCase = null;
        if(PiP!=null){
            PiPCase = branchMap.get(PiP.getNumber());
        }

        if(resampleInfectionTimes) {
            // what happens on i's branch

            if (iCase != iPCase) {
                iCase.setInfectionBranchPosition(MathUtils.nextDouble());
            }

            // what happens between PiP and CiP
            if (PiPCase == null || CiPCase != PiPCase) {
                CiPCase.setInfectionBranchPosition(MathUtils.nextDouble());
            }

        }

        // 3. if the move is down
        if (delta > 0) {


            // 3.1 if the topology will change
            if (PiP != null && tree.getNodeHeight(PiP) < newHeight) {
                // find new parent
                NodeRef newParent = PiP;
                NodeRef newChild = iP;
                while (tree.getNodeHeight(newParent) < newHeight) {
                    newChild = newParent;
                    newParent = tree.getParent(newParent);
                    if (newParent == null) break;
                }


                // if the parent has slid out of its partition
                if(branchMap.get(newChild.getNumber())!=branchMap.get(iP.getNumber())){
//                    throw new OperatorFailedException("invalid slide");
                    return Double.NEGATIVE_INFINITY;
                }

                // if iP is now the earliest node in its subtree


                tree.beginTreeEdit();

                // 3.1.1 if creating a new root
                if (tree.isRoot(newChild)) {
                    tree.removeChild(iP, CiP);
                    tree.removeChild(PiP, iP);
                    tree.addChild(iP, newChild);
                    tree.addChild(PiP, CiP);
                    tree.setRoot(iP);
                    //System.err.println("Creating new root!");

                    if (tree.hasNodeTraits()) {
                        // **********************************************
                        // swap traits and rates so that root keeps it trait and rate values
                        // **********************************************

                        tree.swapAllTraits(newChild, iP);

                    }

                    if (tree.hasRates()) {
                        final double rootNodeRate = tree.getNodeRate(newChild);
                        tree.setNodeRate(newChild, tree.getNodeRate(iP));
                        tree.setNodeRate(iP, rootNodeRate);
                    }

                    // **********************************************

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

                tree.endTreeEdit();


                // 3.1.3 count the hypothetical sources of this destination.
                final int possibleSources = intersectingEdges(tree, newChild, oldHeight, branchMap,
                        branchMap.get(iP.getNumber()), null);
                //System.out.println("possible sources = " + possibleSources);

                logq -= Math.log(possibleSources);

            } else {
                // just change the node height

                tree.setNodeHeight(iP, newHeight);
                logq += 0.0;
            }
        }
        // 4 if we are sliding the subtree up.
        else {

            // 4.0 is it a valid move?
            if (tree.getNodeHeight(i) > newHeight) {
                return Double.NEGATIVE_INFINITY;
            }


            // 4.1 will the move change the topology
            if (tree.getNodeHeight(CiP) > newHeight) {

                List<NodeRef> newChildren = new ArrayList<NodeRef>();
                final int possibleDestinations = intersectingEdges(tree, CiP, newHeight, branchMap,
                        branchMap.get(iP.getNumber()), newChildren);

                // if no valid destinations then return a failure
                if (newChildren.size() == 0) {
                    return Double.NEGATIVE_INFINITY;
                }

                // pick a random parent/child destination edge uniformly from options
                final int childIndex = MathUtils.nextInt(newChildren.size());
                NodeRef newChild = newChildren.get(childIndex);
                NodeRef newParent = tree.getParent(newChild);

                // if newChild was the earliest node in its subtree it has a new infection time

                if(resampleInfectionTimes){
                    AbstractCase newChildCase = branchMap.get(newChild.getNumber());
                    if(newChildCase!=iPCase){
                        newChildCase.setInfectionBranchPosition(MathUtils.nextDouble());
                    }

                }


                tree.beginTreeEdit();

                // 4.1.1 if iP was root
                if (tree.isRoot(iP)) {
                    // new root is CiP
                    tree.removeChild(iP, CiP);
                    tree.removeChild(newParent, newChild);
                    tree.addChild(iP, newChild);
                    tree.addChild(newParent, iP);
                    tree.setRoot(CiP);

                    if (tree.hasNodeTraits()) {
                        // **********************************************
                        // swap traits and rates, so that root keeps it trait and rate values
                        // **********************************************

                        tree.swapAllTraits(iP, CiP);

                    }

                    if (tree.hasRates()) {
                        final double rootNodeRate = tree.getNodeRate(iP);
                        tree.setNodeRate(iP, tree.getNodeRate(CiP));
                        tree.setNodeRate(CiP, rootNodeRate);
                    }

                    // **********************************************

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



                logq += Math.log(possibleDestinations);

            } else {

                tree.setNodeHeight(iP, newHeight);
                logq += 0.0;
            }
        }

        if (swapInRandomRate) {
            final NodeRef j = tree.getNode(MathUtils.nextInt(tree.getNodeCount()));
            if (j != i) {
                final double tmp = tree.getNodeRate(i);
                tree.setNodeRate(i, tree.getNodeRate(j));
                tree.setNodeRate(j, tmp);
            }

        }

        if (swapInRandomTrait) {
            final NodeRef j = tree.getNode(MathUtils.nextInt(tree.getNodeCount()));
            if (j != i) {

                tree.swapAllTraits(i, j);

//                final double tmp = tree.getNodeTrait(i, TRAIT);
//                tree.setNodeTrait(i, TRAIT, tree.getNodeTrait(j, TRAIT));
//                tree.setNodeTrait(j, TRAIT, tmp);
            }

        }

//        if (logq == Double.NEGATIVE_INFINITY) throw new OperatorFailedException("invalid slide");
        if (!Double.isInfinite(logq)) {

            if (DEBUG){
                c2cLikelihood.getTreeModel().checkPartitions();
                c2cLikelihood.outputTreeToFile("afterTSSA.nex", false);
            }

            double reverseEligibleNodeCount = getEligibleNodes(tree, branchMap).size();

            logq += Math.log(eligibleNodeCount/reverseEligibleNodeCount);
        }

        return logq;
    }

    private double getDelta() {
        if (!gaussian) {
            return (MathUtils.nextDouble() * size) - (size / 2.0);
        } else {
            return MathUtils.nextGaussian() * size;
        }
    }

    private boolean eligibleForMove(NodeRef node, TreeModel tree, BranchMapModel branchMap){
        // to be eligible for this move, the node's parent and grandparent, or parent and sibling, must be in the
        // same partition (so removing the parent has no effect on the transmission tree)

        return  (!tree.isRoot(node) && ((tree.getParent(tree.getParent(node))!=null
                && branchMap.get(tree.getParent(node).getNumber())
                ==branchMap.get(tree.getParent(tree.getParent(node)).getNumber()))
                || branchMap.get(tree.getParent(node).getNumber())==branchMap.get(getOtherChild(tree,
                tree.getParent(node), node).getNumber())));
    }

    private ArrayList<NodeRef> getEligibleNodes(TreeModel tree, BranchMapModel branchMap){
        ArrayList<NodeRef> out = new ArrayList<NodeRef>();
        for(NodeRef node : tree.getNodes()){
            if(eligibleForMove(node, tree, branchMap)){
                out.add(node);
            }
        }
        return out;
    }

    //intersectingEdges here is modified to count only possible sources for this special case of the operator - i.e.
    //only branches which have one end in the relevant partition

    private int intersectingEdges(Tree tree, NodeRef node, double height, BranchMapModel branchMap,
                                  AbstractCase partition, List<NodeRef> directChildren) {
        final NodeRef parent = tree.getParent(node);
        if (tree.getNodeHeight(parent) < height || branchMap.get(parent.getNumber())!=partition) return 0;
        if (tree.getNodeHeight(node) < height) {
            if (directChildren != null){
                directChildren.add(node);
            }
            return 1;
        }
        int count = 0;
        for (int i = 0; i < tree.getChildCount(node); i++) {
            count += intersectingEdges(tree, tree.getChild(node, i), height, branchMap, partition, directChildren);
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

    public String getPerformanceSuggestion() {
        return "not implemented";
    }

    public String getOperatorName() {
        return TRANSMISSION_SUBTREE_SLIDE_A + " (" + tree.getId() + ")";
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public static final String RESAMPLE_INFECTION_TIMES = "resampleInfectionTimes";

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            boolean swapRates = xo.getAttribute(SWAP_RATES, false);
            boolean swapTraits = xo.getAttribute(SWAP_TRAITS, false);

            CoercionMode mode = CoercionMode.DEFAULT;
            if (xo.hasAttribute(CoercableMCMCOperator.AUTO_OPTIMIZE)) {
                if (xo.getBooleanAttribute(CoercableMCMCOperator.AUTO_OPTIMIZE)) {
                    mode = CoercionMode.COERCION_ON;
                } else {
                    mode = CoercionMode.COERCION_OFF;
                }
            }

            CaseToCaseTreeLikelihood c2cL = (CaseToCaseTreeLikelihood)xo.getChild(CaseToCaseTreeLikelihood.class);
            final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);

            final double targetAcceptance = xo.getAttribute(TARGET_ACCEPTANCE, 0.234);

            final double size = xo.getAttribute("size", 1.0);

            if (Double.isInfinite(size) || size <= 0.0) {
                throw new XMLParseException("size attribute must be positive and not infinite. was " + size +
                        " for tree " + c2cL.getTreeModel().getId() );
            }

            boolean resampleInfectionTimes = false;

            if(xo.hasAttribute(RESAMPLE_INFECTION_TIMES)) {
                resampleInfectionTimes = xo.getBooleanAttribute(RESAMPLE_INFECTION_TIMES);
            }

            final boolean gaussian = xo.getBooleanAttribute("gaussian");
            TransmissionSubtreeSlideA operator = new TransmissionSubtreeSlideA(c2cL, weight, size, gaussian,
                    swapRates, swapTraits, mode, resampleInfectionTimes);
            operator.setTargetAcceptanceProbability(targetAcceptance);

            return operator;
        }

        public String getParserDescription() {
            return "An operator that slides a phylogenetic subtree, preserving the transmission tree topology.";
        }

        public Class getReturnType() {
            return TransmissionSubtreeSlideA.class;
        }

        public String getParserName() {
            return TRANSMISSION_SUBTREE_SLIDE_A;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
                // Make size optional. If not given or equals zero, size is set to half of average tree branch length.
                AttributeRule.newDoubleRule("size", true),
                AttributeRule.newDoubleRule(TARGET_ACCEPTANCE, true),
                AttributeRule.newBooleanRule("gaussian"),
                AttributeRule.newBooleanRule(SWAP_RATES, true),
                AttributeRule.newBooleanRule(SWAP_TRAITS, true),
                AttributeRule.newBooleanRule(CoercableMCMCOperator.AUTO_OPTIMIZE, true),
                AttributeRule.newBooleanRule(RESAMPLE_INFECTION_TIMES, true),
                new ElementRule(CaseToCaseTreeLikelihood.class)
        };
    };

}
