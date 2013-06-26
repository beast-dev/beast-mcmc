package dr.evomodel.epidemiology.casetocase;

import dr.evolution.tree.NodeRef;
import dr.evomodel.operators.AbstractTreeOperator;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.OperatorFailedException;
import dr.math.MathUtils;
import dr.xml.*;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;

/**
 * Implements branch exchange operations that leave the transmission network unchanged. As this already severely
 * restricts the set of eligible pairs of nodes, this is set up as special case of Wide Exchange.
 *
 * @author Matthew Hall
 */

public class TransmissionExchangeOperatorA extends AbstractTreeOperator {

    private final CaseToCaseTransmissionLikelihood c2cLikelihood;
    public static final String TRANSMISSION_EXCHANGE_OPERATOR_A = "transmissionExchangeOperatorA";

    public TransmissionExchangeOperatorA(CaseToCaseTransmissionLikelihood c2cLikelihood, double weight) {
        this.c2cLikelihood = c2cLikelihood;
        setWeight(weight);
    }

    public double doOperation() throws OperatorFailedException {
        TreeModel tree = c2cLikelihood.getTree();

        final int tipCount = tree.getExternalNodeCount();

        assert tree.getExternalNodeCount() == tipCount : "Lost some tips";

        return 0;
    }

    public void exchange() throws OperatorFailedException{
        TreeModel tree = c2cLikelihood.getTree();

        final int nodeCount = tree.getNodeCount();
        final NodeRef root = tree.getRoot();

        NodeRef i = root;
        NodeRef iP = tree.getParent(i);
        Integer[] possibleSwaps = c2cLikelihood.samePainting(iP,false);
        int noPossibleSwaps = 0;

        while(root == i || noPossibleSwaps == 1) {
            i = tree.getNode(MathUtils.nextInt(nodeCount));
            iP = tree.getParent(i);
            possibleSwaps = c2cLikelihood.samePainting(iP,false);
            noPossibleSwaps = possibleSwaps.length;
        }

        NodeRef jP=iP;

        while(jP==iP){
            jP = tree.getNode(possibleSwaps[MathUtils.nextInt(noPossibleSwaps)]);
        }

        NodeRef j = tree.getChild(jP, MathUtils.nextInt(1));

/*
        I tend to think that this may fail quite a lot of the time due to lack of candidates... a version that does
        actually adjust heights might be necessary in the long run. Narrow exchange might be much more likely to
        actually succeed in changing the tree if the paintings allow the tree to be changed in that way; might
        have to investigate which problem is more serious.
*/

        if( (i != jP) && (j != iP)
                && (tree.getNodeHeight(j) < tree.getNodeHeight(iP))
                && (tree.getNodeHeight(i) < tree.getNodeHeight(jP)) ) {
            exchangeNodes(tree, i, j, iP, jP);
            return;
        }

        c2cLikelihood.makeDirty();

        throw new OperatorFailedException("Couldn't find valid wide move on this tree!");

    }



    public String getPerformanceSuggestion() {
        return "Not implemented";
    }

    public String getOperatorName() {
        return "Transmission tree exchange operator type A (" + c2cLikelihood.getTree().getId() +")";
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return TRANSMISSION_EXCHANGE_OPERATOR_A;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            final CaseToCaseTransmissionLikelihood c2cL
                    = (CaseToCaseTransmissionLikelihood) xo.getChild(CaseToCaseTransmissionLikelihood.class);
            if (c2cL.getTree().getExternalNodeCount() <= 2) {
                throw new XMLParseException("Tree with fewer than 3 taxa");
            }
            final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);

            return new TransmissionExchangeOperatorA(c2cL, weight);
        }

        // ************************************************************************
        // AbstractXMLObjectParser implementation
        // ************************************************************************

        public String getParserDescription(){
            return "This element represents a exchange operator, swapping two random subtrees in such a way that the" +
                    "transmission tree is unaffected.";
        }

        public Class getReturnType(){
            return TransmissionExchangeOperatorA.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules;{
            rules = new XMLSyntaxRule[]{
                    AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
                    new ElementRule(CaseToCaseTransmissionLikelihood.class)
            };
        }
    };

}
