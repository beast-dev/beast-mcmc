package dr.evomodel.epidemiology.casetocase;

import dr.evolution.tree.NodeRef;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;
import dr.xml.*;

import java.util.HashSet;

/**
 * This operator first randomly selects an appropriate node. If its assignment to a case has not previously been changed
 * by this procedure, it then finds its first ancestor which has a different assigned case, and assigns that case to the
 * original node and all intervening nodes. If the assignment had previously been changed, it reverts its assignment
 * and the assignment of all descendants with the same assignment to the assignment of a suitable other descendant.
 *
 * Only necessary for the "extended" case, i.e. where the time of transmission can be later than the TMRCA.
 *
 * @author Matthew Hall
 */
public class NodePaintingCreepOperator extends SimpleMCMCOperator {

    public static final String NODE_PAINTING_CREEP_OPERATOR = "nodePaintingCreepOperator";
    private CaseToCaseTransmissionLikelihood c2cLikelihood;

    public NodePaintingCreepOperator(CaseToCaseTransmissionLikelihood c2cLikelihood, double weight){
        this.c2cLikelihood = c2cLikelihood;
        setWeight(weight);
    }

    public String getOperatorName(){
        return NODE_PAINTING_CREEP_OPERATOR;
    }

    public String getPerformanceSuggestion(){
        return "Not implemented";
    }

    public double doOperation() throws OperatorFailedException {
        int oldUnlockedNodes = 0;
        TreeModel tree = c2cLikelihood.getTree();
        for(int i=0; i<tree.getInternalNodeCount(); i++){
            if(c2cLikelihood.getCreepLocks()[i]){
                oldUnlockedNodes++;
            }
        }
        int internalNodeCount = tree.getInternalNodeCount();
        int nodeToSwitch = MathUtils.nextInt(internalNodeCount);
        // Cannot apply this operator if the node is locked to it
        while(c2cLikelihood.getCreepLocks()[nodeToSwitch]){
            nodeToSwitch = MathUtils.nextInt(internalNodeCount);
        }
        boolean tipLinked = c2cLikelihood.tipLinked(tree.getInternalNode(nodeToSwitch));
        NodeRef node = tree.getInternalNode(nodeToSwitch);
        adjustTree(tree, node, c2cLikelihood.getBranchMap());
        c2cLikelihood.recalculateLocks();
        int newUnlockedNodes = 0;
        for(int i=0; i<tree.getInternalNodeCount(); i++){
            if(c2cLikelihood.getCreepLocks()[i]){
                newUnlockedNodes++;
            }
        }
        if(tipLinked){
            return Math.log(oldUnlockedNodes/(newUnlockedNodes*2));
        } else {
            return Math.log((newUnlockedNodes * 2)/oldUnlockedNodes);
        }
    }

    private void adjustTree(TreeModel tree, NodeRef node, AbstractCase[] map){
        AbstractCase currentCase = map[node.getNumber()];
        if(c2cLikelihood.tipLinked(node)){
            HashSet<Integer> nodesToChange = c2cLikelihood.samePaintingUpTree(node, true);
            NodeRef currentAncestor = node;
            // Find the painting of the first ancestor that doesn't have the same painting
            while(map[currentAncestor.getNumber()]==currentCase){
                currentAncestor = tree.getParent(node);
            }
            // Make the changes
            AbstractCase newCase = map[currentAncestor.getNumber()];
            map[node.getNumber()]=newCase;
            for(int i:nodesToChange){
                map[i]=newCase;
            }
        } else {
            HashSet<Integer> nodesToChange = c2cLikelihood.samePaintingDownTree(node, true);
            nodesToChange.add(node.getNumber());
            NodeRef descendant = node;
            while(c2cLikelihood.countChildrenWithSamePainting(descendant)!=0){
                if(c2cLikelihood.countChildrenWithSamePainting(descendant)>1){
                    throw new RuntimeException("A node that should be creep-locked is not.");
                } else {
                    for(int i=0; i<tree.getChildCount(descendant);i++){
                        if(map[tree.getChild(descendant,i).getNumber()]==currentCase){
                            descendant=tree.getChild(descendant,i);
                        }
                    }
                }
            }
            int choice = MathUtils.nextInt(1);
            AbstractCase replacementPainting = map[tree.getChild(descendant,choice).getNumber()];
            for(Integer i:nodesToChange){
                map[i]=replacementPainting;
            }
        }
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser(){

        public String getParserName(){
            return NODE_PAINTING_CREEP_OPERATOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            CaseToCaseTransmissionLikelihood ftLikelihood =
                    (CaseToCaseTransmissionLikelihood) xo.getChild(CaseToCaseTransmissionLikelihood.class);

            if(!ftLikelihood.isExtended()) {
                throw new XMLParseException("Only extended node paintings use the creep operator.");
            }

            final double weight = xo.getDoubleAttribute("weight");
            return new NodePaintingSwitchOperator(ftLikelihood, weight);
        }

        public String getParserDescription(){
            return "This operator switches the painting of a random eligible internal node from the painting of its " +
                    "first ancestor that has a different painting, or reverses this process";
        }

        public Class getReturnType() {
            return NodePaintingCreepOperator.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newDoubleRule("weight"),
                new ElementRule(CaseToCaseTransmissionLikelihood.class),
        };
    };

}
