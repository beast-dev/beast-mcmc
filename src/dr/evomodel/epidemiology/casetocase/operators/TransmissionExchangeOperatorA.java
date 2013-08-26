package dr.evomodel.epidemiology.casetocase.operators;

import dr.evolution.tree.NodeRef;
import dr.evomodel.epidemiology.casetocase.CaseToCaseTreeLikelihood;
import dr.evomodel.operators.AbstractTreeOperator;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.OperatorFailedException;
import dr.math.MathUtils;
import dr.xml.*;

import java.util.ArrayList;

/**
 * Implements branch exchange operations that leave the transmission network unchanged. As this already severely
 * restricts the set of eligible pairs of nodes, this is set up as special case of Wide Exchange.

 *
 * @author Matthew Hall
 */

public class TransmissionExchangeOperatorA extends AbstractTreeOperator {

    private final CaseToCaseTreeLikelihood c2cLikelihood;
    public static final String TRANSMISSION_EXCHANGE_OPERATOR_A = "transmissionExchangeOperatorA";

    public TransmissionExchangeOperatorA(CaseToCaseTreeLikelihood c2cLikelihood, double weight) {
        this.c2cLikelihood = c2cLikelihood;
        setWeight(weight);
    }

    public double doOperation() throws OperatorFailedException {
        TreeModel tree = c2cLikelihood.getTree();

        double hr = exchange();

        final int tipCount = tree.getExternalNodeCount();

        assert tree.getExternalNodeCount() == tipCount : "Lost some tips";

        return hr;
    }

    public double exchange() throws OperatorFailedException{
        TreeModel tree = c2cLikelihood.getTree();

        final int nodeCount = tree.getNodeCount();
        final NodeRef root = tree.getRoot();

        NodeRef i = root;

        // find a node

        while(root == i){
            i = tree.getNode(MathUtils.nextInt(nodeCount));
        }

        ArrayList<NodeRef> candidates = getPossibleExchanges(tree, i);

        int candidateCount = candidates.size();

        if(candidateCount==0){
            throw new OperatorFailedException("No valid exchanges for this node");
        }

        NodeRef j = candidates.get(MathUtils.nextInt(candidates.size()));

        int jFirstCandidateCount = getPossibleExchanges(tree, j).size();

        double HRDenom = (1/candidateCount) + (1/jFirstCandidateCount);

        NodeRef iP = tree.getParent(i);
        NodeRef jP = tree.getParent(j);

/*
        I tend to think that this may fail quite a lot of the time due to lack of candidates... a version that does
        actually adjust heights might be necessary in the long run. Narrow exchange might be much more likely to
        actually succeed in changing the tree if the paintings allow the tree to be changed in that way; might
        have to investigate which problem is more serious.
*/

        exchangeNodes(tree, i, j, iP, jP);

        ArrayList<NodeRef> reverseCandidatesIfirst = getPossibleExchanges(tree, i);
        ArrayList<NodeRef> reverseCandidatesJfirst = getPossibleExchanges(tree, j);

        double HRNum = (1/reverseCandidatesIfirst.size()) + (1/reverseCandidatesJfirst.size());

        return Math.log(HRNum/HRDenom);
    }

    // get a set of candidates for an exchange of a given node. Does not include the node itself or its sibling, so
    // the check for a failed move will be if this set is of size 0

    public ArrayList<NodeRef> getPossibleExchanges(TreeModel tree, NodeRef node){
        ArrayList<NodeRef> out = new ArrayList<NodeRef>();
        NodeRef parent = tree.getParent(node);
        if(parent==null){
            throw new RuntimeException("Can't exchange the root node");
        }
        Integer[] possibleParentSwaps = c2cLikelihood.samePartition(parent, false);
        for(Integer index: possibleParentSwaps){
            NodeRef newParent = tree.getNode(index);
            if(!tree.isExternal(newParent) && newParent!=parent){
                for(int i=0; i<2; i++){
                    NodeRef candidate = tree.getChild(newParent, i);
                    if(candidate != parent
                            && node != newParent
                            && tree.getNodeHeight(candidate) < tree.getNodeHeight(parent)
                            && tree.getNodeHeight(node) < tree.getNodeHeight(newParent)){
                        if(out.contains(candidate) || candidate==node){
                            throw new RuntimeException("Adding a candidate that already exists in the list or" +
                                    " the node itself");
                        }
                        out.add(candidate);
                    }
                }
            }
        }
        return out;
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
            final CaseToCaseTreeLikelihood c2cL
                    = (CaseToCaseTreeLikelihood) xo.getChild(CaseToCaseTreeLikelihood.class);
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
                    new ElementRule(CaseToCaseTreeLikelihood.class)
            };
        }
    };

}
