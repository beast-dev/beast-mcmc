package dr.evomodel.epidemiology.casetocase.operators;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.epidemiology.casetocase.*;
import dr.inference.operators.GibbsOperator;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;
import dr.xml.*;

import java.util.ArrayList;
import java.util.Arrays;
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

        ArrayList<AbstractCase> caseList = new ArrayList<AbstractCase>(c2cTreeLikelihood.getOutbreak().getCases());

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

        AbstractCase anInfector = c2cTreeLikelihood.getInfector(aCase);

        if(anInfector==null){
            // can't move the root case - there must be one and if all the other outbreak' infection branches are known
            // then this must be it
            return;
        }

        PartitionedTreeModel tree = c2cTreeLikelihood.getTreeModel();

        BranchMapModel originalBranchMap = c2cTreeLikelihood.getBranchMap();

        NodeRef tip1 = c2cTreeLikelihood.getTip(aCase);
        NodeRef tip2 = c2cTreeLikelihood.getTip(anInfector);

        NodeRef mrca = Tree.Utils.getCommonAncestor(tree, tip1, tip2);

        ArrayList<NodeRef> leftBridge = new ArrayList<NodeRef>();
        NodeRef currentNode = tip1;

        while(currentNode!=mrca){
            leftBridge.add(currentNode);
            currentNode = c2cTreeLikelihood.getTreeModel().getParent(currentNode);
        }

        ArrayList<NodeRef> rightBridge = new ArrayList<NodeRef>();
        currentNode = tip2;

        while(currentNode!=mrca){
            rightBridge.add(currentNode);
            currentNode = c2cTreeLikelihood.getTreeModel().getParent(currentNode);
        }

        double[] logProbabilities = new double[leftBridge.size()+rightBridge.size()];

        AbstractCase[][] branchMaps = new AbstractCase[leftBridge.size()+rightBridge.size()][];

        // left bridge

        HashSet<Integer> nodesToChange = c2cTreeLikelihood.samePartitionDownTree(tip1, false);

        AbstractCase[] tempBranchMap = Arrays.copyOf(originalBranchMap.getArrayCopy(), tree.getNodeCount());

        for(Integer number : nodesToChange){
            if(!tree.isExternal(tree.getNode(number))){
                tempBranchMap[number] = anInfector;
            }
        }

        // at this point only the tip is in its partition. Step-by-step, add the left bridge.

        for(int i=0; i<leftBridge.size(); i++){
            NodeRef node = leftBridge.get(i);

            if(i>0){
                tempBranchMap[node.getNumber()]=aCase;
                HashSet<Integer> nodesToChangeUp = c2cTreeLikelihood.samePartitionUpTree(node, false);
                for(Integer number : nodesToChangeUp){
                    tempBranchMap[number]=aCase;
                }
            }

            branchMaps[i] = Arrays.copyOf(tempBranchMap, tempBranchMap.length);

            logProbabilities[i] = c2cTransLikelihood.calculateTempLogLikelihood(tempBranchMap);

        }

        // right bridge

        nodesToChange = c2cTreeLikelihood.samePartitionDownTree(tip2, false);

        tempBranchMap = Arrays.copyOf(originalBranchMap.getArrayCopy(), tree.getNodeCount());

        for(Integer number : nodesToChange){
            if(!tree.isExternal(tree.getNode(number))){
                tempBranchMap[number] = aCase;
            }
        }

        for(int i=0; i<rightBridge.size(); i++){
            NodeRef node = rightBridge.get(i);

            if(i>0){
                tempBranchMap[node.getNumber()]=anInfector;
                HashSet<Integer> nodesToChangeUp = c2cTreeLikelihood.samePartitionUpTree(node, false);
                for(Integer number : nodesToChangeUp){
                    tempBranchMap[number]=anInfector;
                }
            }

            branchMaps[branchMaps.length-1-i] = Arrays.copyOf(tempBranchMap, tempBranchMap.length);

            logProbabilities[branchMaps.length-1-i] = c2cTransLikelihood.calculateTempLogLikelihood(tempBranchMap);

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
