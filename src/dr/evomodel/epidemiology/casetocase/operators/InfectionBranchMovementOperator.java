package dr.evomodel.epidemiology.casetocase.operators;

import dr.evolution.tree.NodeRef;
import dr.evomodel.epidemiology.casetocase.AbstractCase;
import dr.evomodel.epidemiology.casetocase.BranchMapModel;
import dr.evomodel.epidemiology.casetocase.CaseToCaseTreeLikelihood;
import dr.evomodel.epidemiology.casetocase.PartitionedTreeModel;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;
import dr.xml.*;

import java.util.HashSet;

/**
 * This operator finds a branch that corresponds to a transmission event, and moves that event up one branch or down
 * one branch
 *
 * @author Matthew Hall
 */
public class InfectionBranchMovementOperator extends SimpleMCMCOperator{

    public static final String INFECTION_BRANCH_MOVEMENT_OPERATOR = "infectionBranchMovementOperator";
    private CaseToCaseTreeLikelihood c2cLikelihood;
    static final boolean DEBUG = false;

    public InfectionBranchMovementOperator(CaseToCaseTreeLikelihood c2cLikelihood, double weight){
        this.c2cLikelihood = c2cLikelihood;
        setWeight(weight);
    }

    public String getOperatorName(){
        return INFECTION_BRANCH_MOVEMENT_OPERATOR;
    }

    /*  Switch the partition of a randomly selected internal node from the painting of one of its children to the
    * painting of the other, and adjust the rest of the tree to ensure the result still obeys partition rules.*/

    public double doOperation(){

        if(DEBUG){
            c2cLikelihood.debugOutputTree("before.nex", false);
        }

        PartitionedTreeModel tree = c2cLikelihood.getTreeModel();
        BranchMapModel branchMap = c2cLikelihood.getBranchMap();
        int externalNodeCount = tree.getExternalNodeCount();
        // find a case whose infection event we are going to move about
        int nodeToSwitch = MathUtils.nextInt(externalNodeCount);
        // if the infection event is the seed of the epidemic, we need to try again
        while(branchMap.get(tree.getRoot().getNumber())==branchMap.get(tree.getExternalNode(nodeToSwitch).getNumber())){
            nodeToSwitch = MathUtils.nextInt(externalNodeCount);
        }
        // find the child node of the transmission branch
        NodeRef node = tree.getExternalNode(nodeToSwitch);
        while(branchMap.get(node.getNumber())==branchMap.get(tree.getParent(node).getNumber())){
            node = tree.getParent(node);
        }
        double hr = adjustTree(tree, node, branchMap, true);

        if(DEBUG){
            c2cLikelihood.debugOutputTree("after.nex", false);
        }

        return hr;
    }


    private double adjustTree(PartitionedTreeModel tree, NodeRef node, BranchMapModel map, boolean extended){
        // are we going up or down? If we're not extended then all moves are down. External nodes have to move down.
        double out;
        if(!extended || tree.isExternal(node) || MathUtils.nextBoolean()){
            out = moveDown(tree, node, map, extended);
        } else {
            out = moveUp(tree, node, map);
        }
        if(DEBUG){
            c2cLikelihood.checkPartitions();
        }
        return out;
    }

    private double moveDown(PartitionedTreeModel tree, NodeRef node, BranchMapModel map, boolean extended){

        AbstractCase[] newMap = map.getArrayCopy();

        NodeRef parent = tree.getParent(node);

        assert map.get(parent.getNumber())==map.get(node.getNumber()) : "Partition problem";
        if(!extended || c2cLikelihood.tipLinked(parent)){
            NodeRef grandparent = tree.getParent(parent);
            if(grandparent!=null && map.get(grandparent.getNumber())==map.get(parent.getNumber())){

                for(Integer ancestor: c2cLikelihood.samePartitionDownTree(parent, true)){
                    newMap[ancestor]=map.get(node.getNumber());
                }
                newMap[grandparent.getNumber()]=map.get(node.getNumber());
            }
        } else {
            NodeRef sibling = node;
            for(int i=0; i<tree.getChildCount(parent); i++){
                if(tree.getChild(parent, i)!=node){
                    sibling = tree.getChild(parent, i);
                }
            }
            if(map.get(sibling.getNumber())==map.get(parent.getNumber())){
                for(Integer descendant: c2cLikelihood.samePartitionUpTree(sibling, true)){
                    newMap[descendant]=map.get(node.getNumber());
                }
                newMap[sibling.getNumber()]=map.get(node.getNumber());
            }
        }
        newMap[parent.getNumber()]=map.get(node.getNumber());
        map.setAll(newMap);

        return tree.isExternal(node) ? Math.log(0.5) : 0;
    }

    private double moveUp(PartitionedTreeModel tree, NodeRef node, BranchMapModel map){

        AbstractCase[] newMap = map.getArrayCopy();

        double out = 0;

        NodeRef parent = tree.getParent(node);

        assert map.get(parent.getNumber())==map.get(node.getNumber()) : "Partition problem";
        // check if either child is not tip-linked (at most one is not, and if so it must have been in the same
        // partition as both the other child and 'node')
        for(int i=0; i<tree.getChildCount(node); i++){
            NodeRef child = tree.getChild(node, i);
            if(!c2cLikelihood.tipLinked(child)){
                assert map.get(child.getNumber())==map.get(node.getNumber()) : "Partition problem";
                for(Integer descendant: c2cLikelihood.samePartitionUpTree(child, true)){
                    newMap[descendant]=map.get(parent.getNumber());
                }
                newMap[child.getNumber()]=map.get(parent.getNumber());
            } else if(tree.isExternal(child) && map.get(child.getNumber())==map.get(node.getNumber())){
                // we're moving a transmission event onto a terminal branch and need to adjust the HR accordingly
                out += Math.log(2);
            }
        }
        newMap[node.getNumber()]=map.get(parent.getNumber());
        map.setAll(newMap);

        return out;
    }

    public String getPerformanceSuggestion(){
        return "Not implemented";
    }

    /* Parser */

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser(){

        public String getParserName(){
            return INFECTION_BRANCH_MOVEMENT_OPERATOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            CaseToCaseTreeLikelihood ftLikelihood =
                    (CaseToCaseTreeLikelihood) xo.getChild(CaseToCaseTreeLikelihood.class);
            final double weight = xo.getDoubleAttribute("weight");
            return new InfectionBranchMovementOperator(ftLikelihood, weight);
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
                AttributeRule.newDoubleRule("weight"),
                new ElementRule(CaseToCaseTreeLikelihood.class),
        };
    };
}