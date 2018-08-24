/*
 * TransmissionWilsonBaldingA.java
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
import dr.evomodel.epidemiology.casetocase.AbstractCase;
import dr.evomodel.epidemiology.casetocase.BranchMapModel;
import dr.evomodel.epidemiology.casetocase.CaseToCaseTreeLikelihood;
import dr.evomodel.epidemiology.casetocase.PartitionedTreeModel;
import dr.evomodel.operators.AbstractTreeOperator;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.MCMCOperator;
import dr.math.MathUtils;
import dr.xml.*;

import java.util.ArrayList;
import java.util.HashSet;

/**
 * Implements the Wilson-Balding branch swapping move if it does not change the transmission tree.
 *
 * @author Matthew Hall
 */

public class TransmissionWilsonBaldingA extends AbstractTreeOperator {

    private final CaseToCaseTreeLikelihood c2cLikelihood;
    public static final String TRANSMISSION_WILSON_BALDING_A = "transmissionWilsonBaldingA";
    private double logq;
    private static final boolean DEBUG = false;
    private final int tipCount;

    private final boolean resampleInfectionTimes;

    public TransmissionWilsonBaldingA(CaseToCaseTreeLikelihood c2cLikelihood, double weight,
                                      boolean resampleInfectionTimes) {
        this.c2cLikelihood = c2cLikelihood;
        setWeight(weight);
        tipCount = c2cLikelihood.getTreeModel().getExternalNodeCount();

        this.resampleInfectionTimes = resampleInfectionTimes;
    }

    public double doOperation() {
        proposeTree();
        if (c2cLikelihood.getTreeModel().getExternalNodeCount() != tipCount) {
            int newCount = c2cLikelihood.getTreeModel().getExternalNodeCount();
            throw new RuntimeException("Lost some tips in modified SPR! (" +
                    tipCount + "-> " + newCount + ")");
        }

        return logq;
    }

    public void proposeTree() {
        PartitionedTreeModel tree = c2cLikelihood.getTreeModel();
        BranchMapModel branchMap = c2cLikelihood.getBranchMap();
        NodeRef i;
        double oldMinAge, newMinAge, newRange, oldRange, newAge, q;
        // choose a random node avoiding root, and nodes that are ineligible for this move because they have nowhere to
        // go

        ArrayList<NodeRef> eligibleNodes = getEligibleNodes(tree, branchMap);

        i = eligibleNodes.get(MathUtils.nextInt(eligibleNodes.size()));

        double eligibleNodeCount = eligibleNodes.size();

        final NodeRef iP = tree.getParent(i);

        Integer[] sameElements = tree.samePartitionElement(iP);

        HashSet<Integer> possibleDestinations = new HashSet<Integer>();
        // we can insert the node above OR BELOW any node in the same partition
        for (Integer sameElement : sameElements) {
            possibleDestinations.add(sameElement);
            if (!tree.isExternal(tree.getNode(sameElement))) {
                possibleDestinations.add(tree.getChild(tree.getNode(sameElement), 0).getNumber());
                possibleDestinations.add(tree.getChild(tree.getNode(sameElement), 1).getNumber());
            }
        }
        Integer[] pd = possibleDestinations.toArray(new Integer[possibleDestinations.size()]);

        NodeRef j = tree.getNode(pd[MathUtils.nextInt(pd.length)]);
        NodeRef jP = tree.getParent(j);

        while ((jP != null && (tree.getNodeHeight(jP) <= tree.getNodeHeight(i))) || (i == j)) {
            j = tree.getNode(pd[MathUtils.nextInt(pd.length)]);
            jP = tree.getParent(j);
        }

        if (iP == tree.getRoot() || j == tree.getRoot()) {
            logq = Double.NEGATIVE_INFINITY;
            return;
        }

        if (jP == iP || j == iP || jP == i){
            logq = Double.NEGATIVE_INFINITY;
            return;
        }

        final NodeRef CiP = getOtherChild(tree, iP, i);
        NodeRef PiP = tree.getParent(iP);

        if(resampleInfectionTimes) {

            AbstractCase iCase = branchMap.get(i.getNumber());
            AbstractCase iPCase = branchMap.get(iP.getNumber());
            AbstractCase CiPCase = branchMap.get(CiP.getNumber());
            AbstractCase PiPCase = null;
            if(PiP!=null){
                PiPCase = branchMap.get(PiP.getNumber());
            }

            // what happens on i's branch

            if (iCase != iPCase) {
                iCase.setInfectionBranchPosition(MathUtils.nextDouble());
            }

            // what happens between PiP and CiP
            if (PiPCase == null || CiPCase != PiPCase) {
                CiPCase.setInfectionBranchPosition(MathUtils.nextDouble());
            }

            // what happens between k and j


            AbstractCase jCase = branchMap.get(j.getNumber());
            AbstractCase kCase = branchMap.get(jP.getNumber());

            if(iPCase != jCase && iPCase != kCase){
                throw new RuntimeException("TWBA misbehaving.");
            }

            jCase.setInfectionBranchPosition(MathUtils.nextDouble());

        }

        newMinAge = Math.max(tree.getNodeHeight(i), tree.getNodeHeight(j));
        newRange = tree.getNodeHeight(jP) - newMinAge;
        newAge = newMinAge + (MathUtils.nextDouble() * newRange);
        oldMinAge = Math.max(tree.getNodeHeight(i), tree.getNodeHeight(CiP));
        oldRange = tree.getNodeHeight(PiP) - oldMinAge;
        q = newRange / Math.abs(oldRange);

        tree.beginTreeEdit();

        if (j == tree.getRoot()) {

            // 1. remove edges <iP, CiP>
            tree.removeChild(iP, CiP);
            tree.removeChild(PiP, iP);

            // 2. add edges <k, iP>, <iP, j>, <PiP, CiP>
            tree.addChild(iP, j);
            tree.addChild(PiP, CiP);

            // iP is the new root
            tree.setRoot(iP);

        } else if (iP == tree.getRoot()) {

            // 1. remove edges <k, j>, <iP, CiP>, <PiP, iP>
            tree.removeChild(jP, j);
            tree.removeChild(iP, CiP);

            // 2. add edges <k, iP>, <iP, j>, <PiP, CiP>
            tree.addChild(iP, j);
            tree.addChild(jP, iP);

            //CiP is the new root
            tree.setRoot(CiP);

        } else {
            // 1. remove edges <k, j>, <iP, CiP>, <PiP, iP>
            tree.removeChild(jP, j);
            tree.removeChild(iP, CiP);
            tree.removeChild(PiP, iP);

            // 2. add edges <k, iP>, <iP, j>, <PiP, CiP>
            tree.addChild(iP, j);
            tree.addChild(jP, iP);
            tree.addChild(PiP, CiP);
        }

        tree.setNodeHeight(iP, newAge);

        tree.endTreeEdit();

        if(DEBUG){
            c2cLikelihood.getTreeModel().checkPartitions();
        }
        logq = Math.log(q);

        double reverseEligibleNodeCount = getEligibleNodes(tree, branchMap).size();

        logq += Math.log(eligibleNodeCount/reverseEligibleNodeCount);

    }

    public String getPerformanceSuggestion() {
        return "Not implemented";
    }

    private boolean eligibleForMove(NodeRef node, TreeModel tree, BranchMapModel branchMap){
        // to be eligible for this move, the node's parent and grandparent, or parent and other child, must be in the
        // same partition (so removing the parent has no effect on the transmission tree)
        return  (!tree.isRoot(node)
                && ((tree.getParent(tree.getParent(node))!=null
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

    public String getOperatorName() {
        return TRANSMISSION_WILSON_BALDING_A + " (" + c2cLikelihood.getTreeModel().getId() +")";
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public static final String RESAMPLE_INFECTION_TIMES = "resampleInfectionTimes";

        public String getParserName() {
            return TRANSMISSION_WILSON_BALDING_A;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            final CaseToCaseTreeLikelihood c2cL
                    = (CaseToCaseTreeLikelihood) xo.getChild(CaseToCaseTreeLikelihood.class);
            final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);

            boolean resampleInfectionTimes = false;

            if(xo.hasAttribute(RESAMPLE_INFECTION_TIMES)) {
                resampleInfectionTimes = xo.getBooleanAttribute(RESAMPLE_INFECTION_TIMES);
            }

            return new TransmissionWilsonBaldingA(c2cL, weight, resampleInfectionTimes);
        }

        // ************************************************************************
        // AbstractXMLObjectParser implementation
        // ************************************************************************

        public String getParserDescription(){
            return "This element represents a Wilson-Balding move operator, such that the transplantation of the " +
                    "subtree does not affect the topology of the transmission tree.";
        }

        public Class getReturnType(){
            return TransmissionWilsonBaldingA.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules;{
            rules = new XMLSyntaxRule[]{
                    AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
                    AttributeRule.newBooleanRule(RESAMPLE_INFECTION_TIMES, true),
                    new ElementRule(CaseToCaseTreeLikelihood.class)
            };
        }
    };
}
