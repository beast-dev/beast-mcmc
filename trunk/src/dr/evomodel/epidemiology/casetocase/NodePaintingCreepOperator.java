package dr.evomodel.epidemiology.casetocase;

import dr.evolution.tree.NodeRef;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;
import dr.xml.*;

/**
 * Created with IntelliJ IDEA.
 * User: mhall
 * Date: 24/01/2013
 * Time: 16:18
 * To change this template use File | Settings | File Templates.
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
        TreeModel tree = c2cLikelihood.getTree();
        int internalNodeCount = tree.getInternalNodeCount();
        int nodeToSwitch = MathUtils.nextInt(internalNodeCount);
        while(c2cLikelihood.getCreepLocks()[nodeToSwitch]){
            nodeToSwitch = MathUtils.nextInt(internalNodeCount);
        }
        NodeRef node = tree.getInternalNode(nodeToSwitch);
        AbstractCase currentCase = c2cLikelihood.getBranchMap()[node.getNumber()];
        return 0;
    }

    private void adjustTree(NodeRef node){
        if(c2cLikelihood.tipLinked(node)){

        } else {

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
