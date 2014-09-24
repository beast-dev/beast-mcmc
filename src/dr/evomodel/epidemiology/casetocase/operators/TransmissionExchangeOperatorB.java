package dr.evomodel.epidemiology.casetocase.operators;

import dr.evolution.tree.NodeRef;
import dr.evomodel.epidemiology.casetocase.AbstractCase;
import dr.evomodel.epidemiology.casetocase.AbstractOutbreak;
import dr.evomodel.epidemiology.casetocase.BranchMapModel;
import dr.evomodel.epidemiology.casetocase.CaseToCaseTreeLikelihood;
import dr.evomodel.operators.AbstractTreeOperator;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.OperatorFailedException;
import dr.math.MathUtils;
import dr.xml.*;

import java.util.ArrayList;

/**
 * Implements branch exchange operations that also exchange entire subtrees of the transmission tree. As this already
 * severely restricts the set of eligible pairs of nodes, this is set up as special case of Wide Exchange.
 *
 * @author Matthew Hall
 */

public class TransmissionExchangeOperatorB extends AbstractTreeOperator {

    private final CaseToCaseTreeLikelihood c2cLikelihood;
    public static final String TRANSMISSION_EXCHANGE_OPERATOR_B = "transmissionExchangeOperatorB";

    private final boolean resampleInfectionTimes;

    public TransmissionExchangeOperatorB(CaseToCaseTreeLikelihood c2cLikelihood, double weight,
                                         boolean resampleInfectionTimes) {
        this.c2cLikelihood = c2cLikelihood;
        setWeight(weight);

        this.resampleInfectionTimes = resampleInfectionTimes;
    }

    public double doOperation() throws OperatorFailedException {

        TreeModel tree = c2cLikelihood.getTreeModel();

        double hr = exchange();

        final int tipCount = tree.getExternalNodeCount();

        assert tree.getExternalNodeCount() == tipCount : "Lost some tips";

        return hr;
    }

    public double exchange() throws OperatorFailedException{
        TreeModel tree = c2cLikelihood.getTreeModel();
        BranchMapModel branchMap = c2cLikelihood.getBranchMap();

        final int nodeCount = tree.getNodeCount();
        final NodeRef root = tree.getRoot();

        NodeRef i = root;
        NodeRef iP = tree.getParent(i);
        boolean partitionsMatch = true;

        // find a node - its parent must be in a different partition (this will not be different following the move)

        while(root == i || partitionsMatch){
            i = tree.getNode(MathUtils.nextInt(nodeCount));
            iP = tree.getParent(i);
            partitionsMatch = i == root || branchMap.get(i.getNumber()) == branchMap.get(iP.getNumber());
        }

        ArrayList<NodeRef> candidates = getPossibleExchanges(tree, i);

        int candidateCount = candidates.size();

        if(candidateCount==0){
            throw new OperatorFailedException("No valid exchanges for this node");
        }

        NodeRef j = candidates.get(MathUtils.nextInt(candidates.size()));

        int jFirstCandidateCount = getPossibleExchanges(tree, j).size();

        double HRDenom = (1/((double)candidateCount)) + (1/((double)jFirstCandidateCount));

        NodeRef jP = tree.getParent(j);


        if(resampleInfectionTimes){
            Parameter branchPostitions = c2cLikelihood.getInfectionTimeBranchPositions();

            AbstractCase iCase = branchMap.get(i.getNumber());
            AbstractCase jCase = branchMap.get(j.getNumber());
            AbstractCase iPCase = branchMap.get(iP.getNumber());
            AbstractCase jPCase = branchMap.get(jP.getNumber());

            branchPostitions.setParameterValue(c2cLikelihood.getOutbreak().getCaseIndex(iCase),
                    MathUtils.nextDouble());
            branchPostitions.setParameterValue(c2cLikelihood.getOutbreak().getCaseIndex(jCase),
                    MathUtils.nextDouble());

        }

/*
        Intuitively it would seem this is a lot more likely to succeed than operator A.
*/

        exchangeNodes(tree, i, j, iP, jP);

        ArrayList<NodeRef> reverseCandidatesIfirst = getPossibleExchanges(tree, i);
        ArrayList<NodeRef> reverseCandidatesJfirst = getPossibleExchanges(tree, j);

        double HRNum = (1/((double)reverseCandidatesIfirst.size())) + (1/((double)reverseCandidatesJfirst.size()));

        return Math.log(HRNum/HRDenom);

    }

    // get a set of candidates for an exchange of a given node. Does not include the node itself or its sibling, so
    // the check for a failed move will be if this set is of size 0

    public ArrayList<NodeRef> getPossibleExchanges(TreeModel tree, NodeRef node){
        BranchMapModel map = c2cLikelihood.getBranchMap();
        ArrayList<NodeRef> out = new ArrayList<NodeRef>();
        NodeRef parent = tree.getParent(node);
        if(parent==null){
            throw new RuntimeException("Can't exchange the root node");
        }
        if(map.get(parent.getNumber())==map.get(node.getNumber())){
            throw new RuntimeException("This node is not exchangeable by this operator");
        }
        for(NodeRef candidate: tree.getNodes()){
            NodeRef newParent = tree.getParent(candidate);
            if(newParent!=parent && newParent!=null){
                if(candidate != parent
                        && node != newParent
                        && tree.getNodeHeight(candidate) < tree.getNodeHeight(parent)
                        && tree.getNodeHeight(node) < tree.getNodeHeight(newParent)
                        && map.get(newParent.getNumber())!=map.get(candidate.getNumber())){
                    if(out.contains(candidate) || candidate==node){
                        throw new RuntimeException("Adding a candidate that already exists in the list or" +
                                " the node itself");
                    }
                    out.add(candidate);
                }

            }
        }
        return out;
    }



    public String getPerformanceSuggestion() {
        return "Not implemented";
    }

    public String getOperatorName() {
        return "Transmission tree exchange operator type B (" + c2cLikelihood.getTreeModel().getId() +")";
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public static final String RESAMPLE_INFECTION_TIMES = "resampleInfectionTimes";

        public String getParserName() {
            return TRANSMISSION_EXCHANGE_OPERATOR_B;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            final CaseToCaseTreeLikelihood c2cL
                    = (CaseToCaseTreeLikelihood) xo.getChild(CaseToCaseTreeLikelihood.class);
            if (c2cL.getTreeModel().getExternalNodeCount() <= 2) {
                throw new XMLParseException("Tree with fewer than 3 taxa");
            }
            final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);

            boolean resampleInfectionTimes = false;

            if(xo.hasAttribute(RESAMPLE_INFECTION_TIMES)) {
                resampleInfectionTimes = xo.getBooleanAttribute(RESAMPLE_INFECTION_TIMES);
            }

            return new TransmissionExchangeOperatorB(c2cL, weight, resampleInfectionTimes);
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
                    AttributeRule.newBooleanRule(RESAMPLE_INFECTION_TIMES, true),
                    new ElementRule(CaseToCaseTreeLikelihood.class),
            };
        }
    };


}
