package dr.evomodel.epidemiology.casetocase.operators;

import dr.evolution.tree.NodeRef;
import dr.evomodel.epidemiology.casetocase.AbstractCase;
import dr.evomodel.epidemiology.casetocase.CaseToCaseTreeLikelihood;
import dr.evomodel.operators.AbstractTreeOperator;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.OperatorFailedException;
import dr.math.MathUtils;
import dr.xml.*;

/**
 * Implements branch exchange operations that also exchange entire subtrees of the transmission tree. As this already
 * severely restricts the set of eligible pairs of nodes, this is set up as special case of Wide Exchange.
 *
 * @author Matthew Hall
 */

public class TransmissionExchangeOperatorB extends AbstractTreeOperator {

    private final CaseToCaseTreeLikelihood c2cLikelihood;
    public static final String TRANSMISSION_EXCHANGE_OPERATOR_B = "transmissionExchangeOperatorB";

    public TransmissionExchangeOperatorB(CaseToCaseTreeLikelihood c2cLikelihood, double weight) {
        this.c2cLikelihood = c2cLikelihood;
        setWeight(weight);
    }

    public double doOperation() throws OperatorFailedException {
        TreeModel tree = c2cLikelihood.getTree();

        final int tipCount = tree.getExternalNodeCount();

        assert tree.getExternalNodeCount() == tipCount :
                "Lost some tips";

        return 0;
    }

    public void exchange() throws OperatorFailedException{
        TreeModel tree = c2cLikelihood.getTree();
        AbstractCase[] branchMap = c2cLikelihood.getBranchMap();

        final int nodeCount = tree.getNodeCount();
        final NodeRef root = tree.getRoot();

        NodeRef i = root;
        NodeRef iP = tree.getParent(i);
        boolean paintingsDoNotMatch = false;

        while(root == i || !paintingsDoNotMatch) {
            i = tree.getNode(MathUtils.nextInt(nodeCount));
            iP = tree.getParent(i);
            if(root!=i){
                paintingsDoNotMatch = branchMap[i.getNumber()]!=branchMap[iP.getNumber()];
            }
        }

        NodeRef j = root;
        NodeRef jP = tree.getParent(j);
        paintingsDoNotMatch = false;

        while(root == j || i==j || iP==jP || !paintingsDoNotMatch) {
            j = tree.getNode(MathUtils.nextInt(nodeCount));
            jP = tree.getParent(j);
            if(root!=j){
                paintingsDoNotMatch = branchMap[j.getNumber()]!=branchMap[jP.getNumber()];
            }
        }

/*
        Intuitively it would seem this is a lot more likely to succeed than operator A.
*/

        if((i != jP) && (j != iP)
                && (tree.getNodeHeight(j) < tree.getNodeHeight(iP))
                && (tree.getNodeHeight(i) < tree.getNodeHeight(jP)) ) {
            exchangeNodes(tree, i, j, iP, jP);
            return;
        }

        c2cLikelihood.makeDirty();

        throw new OperatorFailedException("Couldn't find valid move on this tree!");

    }



    public String getPerformanceSuggestion() {
        return "Not implemented";
    }

    public String getOperatorName() {
        return "Transmission tree exchange operator type B (" + c2cLikelihood.getTree().getId() +")";
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return TRANSMISSION_EXCHANGE_OPERATOR_B;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            final CaseToCaseTreeLikelihood c2cL
                    = (CaseToCaseTreeLikelihood) xo.getChild(CaseToCaseTreeLikelihood.class);
            if (c2cL.getTree().getExternalNodeCount() <= 2) {
                throw new XMLParseException("Tree with fewer than 3 taxa");
            }
            final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);

            return new TransmissionExchangeOperatorB(c2cL, weight);
        }

        // ************************************************************************
        // AbstractXMLObjectParser implementation
        // ************************************************************************

        public String getParserDescription(){
            return "This element represents a exchange operator, swapping two random subtrees in such a way that " +
                    "subtrees of the transmission tree are also exchanged.";
        }

        public Class getReturnType(){
            return TransmissionExchangeOperatorB.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules;{
            rules = new XMLSyntaxRule[]{
                    AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
                    new ElementRule(CaseToCaseTreeLikelihood.class)
            };
        }
    };


}
