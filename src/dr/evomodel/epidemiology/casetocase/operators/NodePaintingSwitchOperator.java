package dr.evomodel.epidemiology.casetocase.operators;

import dr.evolution.tree.NodeRef;
import dr.evomodel.epidemiology.casetocase.AbstractCase;
import dr.evomodel.epidemiology.casetocase.CaseToCaseTreeLikelihood;
import dr.evomodel.tree.TreeModel;
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
public class NodePaintingSwitchOperator extends SimpleMCMCOperator{

    public static final String NODE_PAINTING_SWITCH_OPERATOR = "nodePaintingSwitchOperator";
    private CaseToCaseTreeLikelihood c2cLikelihood;
    boolean debug = false;

    public NodePaintingSwitchOperator(CaseToCaseTreeLikelihood c2cLikelihood, double weight){
        this.c2cLikelihood = c2cLikelihood;
        setWeight(weight);
    }

    public String getOperatorName(){
        return NODE_PAINTING_SWITCH_OPERATOR;
    }

    /*  Switch the painting of a randomly selected internal node from the painting of one of its children to the
    * painting of the other, and adjust the rest of the tree to ensure the result is still a painting.*/

    public double doOperation(){
        TreeModel tree = c2cLikelihood.getTree();
        AbstractCase[] branchMap = c2cLikelihood.getBranchMap();
        int externalNodeCount = tree.getExternalNodeCount();
        // find a case whose infection event we are going to move about
        int nodeToSwitch = MathUtils.nextInt(externalNodeCount);
        // if the infection event is the seed of the epidemic, we need to try again
        while(branchMap[tree.getRoot().getNumber()]==branchMap[tree.getExternalNode(nodeToSwitch).getNumber()]){
            nodeToSwitch = MathUtils.nextInt(externalNodeCount);
        }
        // find the child node of the transmission branch
        NodeRef node = tree.getExternalNode(nodeToSwitch);
        while(branchMap[node.getNumber()]==branchMap[tree.getParent(node).getNumber()]){
            node = tree.getParent(node);
        }
        double hr = adjustTree(tree, node, branchMap, c2cLikelihood.getRecalculationArray(), c2cLikelihood.isExtended());
        c2cLikelihood.makeDirty(false);
        return hr;
    }


    private double adjustTree(TreeModel tree, NodeRef node, AbstractCase[] map, boolean[] recalcArray, boolean extended){
        // are we going up or down? If we're not extended then all moves are down. External nodes have to move down.
        double out;
        if(!extended || tree.isExternal(node) || MathUtils.nextBoolean()){
            out = moveDown(tree, node, map, extended);
        } else {
            out = moveUp(tree, node, map);
        }
        if(debug){
            c2cLikelihood.checkPartitions();
        }
        return out;
    }

    private double moveDown(TreeModel tree, NodeRef node, AbstractCase[] map, boolean extended){
        NodeRef parent = tree.getParent(node);
        assert map[parent.getNumber()]==map[node.getNumber()] : "Partition problem";
        if(!extended || c2cLikelihood.tipLinked(parent)){
            NodeRef grandparent = tree.getParent(parent);
            if(grandparent!=null && map[grandparent.getNumber()]==map[parent.getNumber()]){
                for(Integer ancestor: c2cLikelihood.samePartitionDownTree(parent, true)){
                    map[ancestor]=map[node.getNumber()];
                }
                map[grandparent.getNumber()]=map[node.getNumber()];
            }
        } else {
            NodeRef sibling = node;
            for(int i=0; i<tree.getChildCount(parent); i++){
                if(tree.getChild(parent, i)!=node){
                    sibling = tree.getChild(parent, i);
                }
            }
            if(map[sibling.getNumber()]==map[parent.getNumber()]){
                for(Integer descendant: c2cLikelihood.samePartitionUpTree(sibling, true)){
                    map[descendant]=map[node.getNumber()];
                }
                map[sibling.getNumber()]=map[node.getNumber()];
            }
        }
        map[parent.getNumber()]=map[node.getNumber()];
        c2cLikelihood.flagForDescendantRecalculation(tree, node);
        return tree.isExternal(node) ? Math.log(0.5) : 0;
    }

    private double moveUp(TreeModel tree, NodeRef node, AbstractCase[] map){
        double out = 0;
        NodeRef parent = tree.getParent(node);
        assert map[parent.getNumber()]==map[node.getNumber()] : "Partition problem";
        // check if either child is not tip-linked (at most one is not, and if so it must have been in the same
        // partition as both the other child and 'node')
        for(int i=0; i<tree.getChildCount(node); i++){
            NodeRef child = tree.getChild(node, i);
            if(!c2cLikelihood.tipLinked(child)){
                assert map[child.getNumber()]==map[node.getNumber()] : "Partition problem";
                for(Integer descendant: c2cLikelihood.samePartitionUpTree(child, true)){
                    map[descendant]=map[parent.getNumber()];
                }
                map[child.getNumber()]=map[parent.getNumber()];
            } else if(tree.isExternal(child) && map[child.getNumber()]==map[node.getNumber()]){
                // we're moving a transmission event onto a terminal branch and need to adjust the HR accordingly
                out += Math.log(2);
            }
        }
        map[node.getNumber()]=map[parent.getNumber()];
        return out;
    }

    public String getPerformanceSuggestion(){
        return "Not implemented";
    }

    /* Parser */

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser(){

        public String getParserName(){
            return NODE_PAINTING_SWITCH_OPERATOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            CaseToCaseTreeLikelihood ftLikelihood =
                    (CaseToCaseTreeLikelihood) xo.getChild(CaseToCaseTreeLikelihood.class);
            final double weight = xo.getDoubleAttribute("weight");
            return new NodePaintingSwitchOperator(ftLikelihood, weight);
        }

        public String getParserDescription(){
            return "This operator switches the painting of a random eligible internal node from the painting of one of " +
                    "its children to the painting of the other";
        }

        public Class getReturnType() {
            return NodePaintingSwitchOperator.class;
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