/*
 * InfectionBranchMovementOperator.java
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
import dr.evomodel.epidemiology.casetocase.*;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;
import dr.xml.*;

/**
 * This operator finds a branch that corresponds to a transmission event, and moves that event up one branch or down
 * one branch
 *
 * @author Matthew Hall
 */
public class InfectionBranchMovementOperator extends SimpleMCMCOperator{

    public static final String INFECTION_BRANCH_MOVEMENT_OPERATOR = "infectionBranchMovementOperator";
    private CaseToCaseTreeLikelihood c2cLikelihood;

    private final boolean resampleInfectionTimes;

    private static final boolean DEBUG = false;

    public InfectionBranchMovementOperator(CaseToCaseTreeLikelihood c2cLikelihood, double weight,
                                           boolean resampleInfectionTimes){
        this.c2cLikelihood = c2cLikelihood;
        setWeight(weight);

        this.resampleInfectionTimes = resampleInfectionTimes;
    }

    public String getOperatorName(){
        return INFECTION_BRANCH_MOVEMENT_OPERATOR;
    }

    /*  Switch the partition of a randomly selected internal node from the painting of one of its children to the
    * painting of the other, and adjust the rest of the tree to ensure the result still obeys partition rules.*/

    public double doOperation() {

        PartitionedTreeModel tree = c2cLikelihood.getTreeModel();
        BranchMapModel branchMap = c2cLikelihood.getBranchMap();

        // find a case whose infection event we are going to move about
        int caseIndexToAdjust = MathUtils.nextInt(c2cLikelihood.getOutbreak().size());

        AbstractCase aCase = c2cLikelihood.getOutbreak().getCase(caseIndexToAdjust);

        // if the infection event is the seed of the epidemic, we need to try again
        while(branchMap.get(tree.getRoot().getNumber()) == aCase || !aCase.wasEverInfected()){
            caseIndexToAdjust = MathUtils.nextInt(c2cLikelihood.getOutbreak().size());
            aCase = c2cLikelihood.getOutbreak().getCase(caseIndexToAdjust);
        }

        // find the child node of the transmission branch
        NodeRef node = tree.getEarliestNodeInElement(aCase);

        double hr = adjustTree(tree, node);

        return hr;
    }


    private double adjustTree(PartitionedTreeModel tree, NodeRef node) {
        double out;

        BranchMapModel map = tree.getBranchMap();

        AbstractCase infectedCase = map.get(node.getNumber());
        AbstractCase infectorCase = map.get(tree.getParent(node).getNumber());

        NodeRef infectedMRCA = tree.caseMRCA(infectedCase);

        boolean downIsPossible = node != infectedMRCA;
        boolean upIsPossible = !(tree.isRootBlockedBy(infectedCase, infectorCase)
                && tree.isAncestral(tree.getParent(node)));

        if(upIsPossible && downIsPossible){
            out = MathUtils.nextBoolean() ? moveUp(tree, node) : moveDown(tree, node);
        } else if(upIsPossible){
            out = moveUp(tree, node);
        } else if(downIsPossible){
            out = moveDown(tree, node);
        } else {
//            throw new OperatorFailedException("Chosen infection event cannot be adjusted in this tree");
            return Double.NEGATIVE_INFINITY;
        }

        if(DEBUG){
            c2cLikelihood.getTreeModel().checkPartitions();
        }
        return out;
    }

    private double moveUp(PartitionedTreeModel tree, NodeRef node){
        BranchMapModel map = tree.getBranchMap();

        AbstractCase infectedCase = map.get(node.getNumber());

        AbstractCase[] newMap = map.getArrayCopy();

        NodeRef parent = tree.getParent(node);

        double hr = 0;

        NodeRef sibling = node;
        for(int i=0; i<tree.getChildCount(parent); i++){
            if(tree.getChild(parent, i)!=node){
                sibling = tree.getChild(parent, i);
            }
        }

        AbstractCase infectorCase = map.get(parent.getNumber());

        NodeRef infectedMRCA = tree.caseMRCA(infectedCase);
        NodeRef infectorMRCA = tree.caseMRCA(infectorCase);

        if(c2cLikelihood.getTreeModel().isAncestral(parent)){

            if(resampleInfectionTimes){
                infectorCase.setInfectionBranchPosition(MathUtils.nextDouble());
            }

            NodeRef grandparent = tree.getParent(parent);
            if(grandparent!=null && map.get(grandparent.getNumber())==map.get(parent.getNumber())){
                for(Integer ancestor: c2cLikelihood.getTreeModel().samePartitionElementUpTree(parent)){
                    newMap[ancestor] = map.get(node.getNumber());
                }
                newMap[grandparent.getNumber()]=map.get(node.getNumber());
            }

            hr += node == infectedMRCA ? Math.log(0.5) : 0;

        } else {
            if(map.get(sibling.getNumber())==map.get(parent.getNumber())){
                for(Integer descendant: c2cLikelihood.getTreeModel().samePartitionElementDownTree(sibling)){
                    newMap[descendant]=map.get(node.getNumber());
                }
                newMap[sibling.getNumber()]=map.get(node.getNumber());
            }

            hr += node == infectedMRCA ? Math.log(0.5) : 0;
        }
        newMap[parent.getNumber()]=map.get(node.getNumber());
        map.setAll(newMap, false);

        //HR adjustments for reverse moves
        if(c2cLikelihood.getTreeModel().isAncestral(parent)){
            hr += sibling == infectorMRCA ? Math.log(2) : 0;
        } else {
            NodeRef grandparent = tree.getParent(parent);

            hr += tree.isRootBlockedBy(infectedCase, infectorCase)
                    && tree.isAncestral(grandparent) ? Math.log(2) : 0;

        }


        if(resampleInfectionTimes){
            infectedCase.setInfectionBranchPosition(MathUtils.nextDouble());
        }

        return hr;
    }

    private double moveDown(PartitionedTreeModel tree, NodeRef node){
        BranchMapModel map = tree.getBranchMap();

        AbstractCase infectedCase = map.get(node.getNumber());
        AbstractCase infectorCase = map.get(tree.getParent(node).getNumber());

        AbstractCase[] newMap = map.getArrayCopy();

        double out = 0;

        NodeRef parent = tree.getParent(node);

        NodeRef infectedMRCA = tree.caseMRCA(infectedCase);

        // check if either child is not ancestral (at most one is not, and if so it must have been in the same
        // partition as both the other child and 'node')
        for(int i=0; i<tree.getChildCount(node); i++){
            NodeRef child = tree.getChild(node, i);
            if(!c2cLikelihood.getTreeModel().isAncestral(child)){
                for(Integer descendant: c2cLikelihood.getTreeModel().samePartitionElementDownTree(child)){
                    newMap[descendant]=map.get(parent.getNumber());
                }
                newMap[child.getNumber()]=map.get(parent.getNumber());
            } else if(child == infectedMRCA && map.get(child.getNumber())==map.get(node.getNumber())){
                // we're moving a transmission event as far down as it can go and need to adjust the HR accordingly
                out += Math.log(2);
            }
        }

        //if you couldn't move it any further up
        out += tree.isRootBlockedBy(infectedCase, infectorCase)
                && tree.isAncestral(parent) ? Math.log(0.5) : 0;

        if(resampleInfectionTimes){
            infectedCase.setInfectionBranchPosition(MathUtils.nextDouble());
        }

        newMap[node.getNumber()]=map.get(parent.getNumber());
        map.setAll(newMap, false);

        return out;
    }

    public String getPerformanceSuggestion(){
        return "Not implemented";
    }

    /* Parser */

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser(){

        public static final String RESAMPLE_INFECTION_TIMES = "resampleInfectionTimes";

        public String getParserName(){
            return INFECTION_BRANCH_MOVEMENT_OPERATOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            CaseToCaseTreeLikelihood ftLikelihood =
                    (CaseToCaseTreeLikelihood) xo.getChild(CaseToCaseTreeLikelihood.class);
            final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);

            boolean resampleInfectionTimes = false;

            if(xo.hasAttribute(RESAMPLE_INFECTION_TIMES)) {
                resampleInfectionTimes = xo.getBooleanAttribute(RESAMPLE_INFECTION_TIMES);
            }

            return new InfectionBranchMovementOperator(ftLikelihood, weight, resampleInfectionTimes);
        }

        public String getParserDescription(){
            return "This operator switches the painting of a random eligible internal node from the painting of one " +
                    "of its children to the painting of the other";
        }

        public Class getReturnType() {
            return InfectionBranchMovementOperator.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
                AttributeRule.newBooleanRule(RESAMPLE_INFECTION_TIMES, true),
                new ElementRule(CaseToCaseTreeLikelihood.class),
        };
    };
}