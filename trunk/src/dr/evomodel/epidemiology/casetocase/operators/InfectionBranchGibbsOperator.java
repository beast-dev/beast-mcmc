package dr.evomodel.epidemiology.casetocase.operators;

import dr.evolution.tree.NodeRef;
import dr.evomodel.epidemiology.casetocase.*;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.GibbsOperator;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;
import dr.xml.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;

/**
 * Goes through every case in turn, except the root case whose infection branch is fixed given all the others, and
 * calculates the probability of each legal branch placement for the infection (if there is more than one). Chooses a
 * new one with proportional probability.
 *
 * @author Matthew Hall
 */
public class InfectionBranchGibbsOperator extends SimpleMCMCOperator implements GibbsOperator {

    public static final String INFECTION_BRANCH_GIBBS_OPERATOR = "infectionBranchGibbsOperator";
    private CaseToCaseTransmissionLikelihood c2cTransLikelihood;
    private CaseToCaseTreeLikelihood c2cTreeLikelihood;
    private static final boolean DEBUG = false;

    public InfectionBranchGibbsOperator(CaseToCaseTransmissionLikelihood c2cTransLikelihood, double weight){
        this.c2cTransLikelihood = c2cTransLikelihood;
        c2cTreeLikelihood = c2cTransLikelihood.getTreeLikelihood();
        setWeight(weight);
    }

    public int getStepCount() {
        return c2cTreeLikelihood.getOutbreak().size()-1;
    }

    public String getPerformanceSuggestion() {
        return null;
    }

    public String getOperatorName() {
        return INFECTION_BRANCH_GIBBS_OPERATOR;
    }

    public double doOperation() throws OperatorFailedException {

        ArrayList<AbstractCase> caseList = c2cTreeLikelihood.getOutbreak().getCases();

        int[] shuffledOrder = MathUtils.shuffled(caseList.size());

        for(int i : shuffledOrder){

            if(DEBUG){
                String caseName = caseList.get(i).getName();
                c2cTreeLikelihood.debugOutputTree(caseName + "_before.nex", false);
            }

            pickBranch(caseList.get(i));
        }

        if(DEBUG){
            c2cTreeLikelihood.checkPartitions();
        }

        return 0;
    }

    private void pickBranch(AbstractCase aCase){

        if(c2cTreeLikelihood.getInfector(aCase)==null){
            // can't move the root case - there must be one and if all the other cases' infection branches are known
            // then this must be it
            return;
        }

        PartitionedTreeModel tree = c2cTreeLikelihood.getTreeModel();

        BranchMapModel originalBranchMap = c2cTreeLikelihood.getBranchMap();

        NodeRef tip = c2cTreeLikelihood.getTip(aCase);

        ArrayList<NodeRef> positionsPlease = new ArrayList<NodeRef>();
        NodeRef currentNode = tip;

        // the parent is a candidate branch unless the current branch map has a different case there and it is an
        // ancestor of the tip corresponding to that case

        while(originalBranchMap.get(currentNode.getNumber())==aCase || !c2cTreeLikelihood.tipLinked(currentNode)){
            positionsPlease.add(currentNode);
            currentNode = c2cTreeLikelihood.getTreeModel().getParent(currentNode);
        }

        AbstractCase infector = originalBranchMap.get(currentNode.getNumber());

        if(positionsPlease.size()==1){
            //nowhere to go
            return;
        }

        double[] logProbabilities = new double[positionsPlease.size()];

        HashSet<Integer> nodesToChange = c2cTreeLikelihood.samePartitionDownTree(tip, false);

        AbstractCase[][] branchMaps = new AbstractCase[positionsPlease.size()][];

        AbstractCase[] tempBranchMap = Arrays.copyOf(originalBranchMap.getArrayCopy(), tree.getNodeCount());

        for(Integer number : nodesToChange){
            tempBranchMap[number] = infector;
        }

        branchMaps[0] = tempBranchMap;
        logProbabilities[0] = c2cTransLikelihood.calculateTempLogLikelihood(tempBranchMap);

        // at this point only the tip is in its partition

        for(int i=1; i<positionsPlease.size(); i++){
            NodeRef node = positionsPlease.get(i);
            AbstractCase[] newBranchMap = Arrays.copyOf(branchMaps[i-1], tempBranchMap.length);
            newBranchMap[node.getNumber()]=aCase;
            HashSet<Integer> nodesToChangeUp = c2cTreeLikelihood.samePartitionUpTree(node, false);
            for(Integer number : nodesToChangeUp){
                if(!tree.isExternal(tree.getNode(number))){
                    newBranchMap[number]=aCase;
                }
            }
            branchMaps[i] = newBranchMap;
            logProbabilities[i] = c2cTransLikelihood.calculateTempLogLikelihood(newBranchMap);
        }

        // this prevents underflow

        int choice = MathUtils.randomChoiceLogPDF(logProbabilities);

        originalBranchMap.setAll(branchMaps[choice], false);
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser(){

        public String getParserName(){
            return INFECTION_BRANCH_GIBBS_OPERATOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            CaseToCaseTransmissionLikelihood ftLikelihood =
                    (CaseToCaseTransmissionLikelihood) xo.getChild(CaseToCaseTransmissionLikelihood.class);
            final double weight = xo.getDoubleAttribute("weight");
            return new InfectionBranchGibbsOperator(ftLikelihood, weight);
        }

        public String getParserDescription(){
            return "A Gibbs sampler on the branches that correspond to the infection of each case";
        }

        public Class getReturnType() {
            return InfectionBranchGibbsOperator.class;
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
