package dr.evomodel.epidemiology.casetocase;

import dr.app.tools.NexusExporter;
import dr.evolution.coalescent.Coalescent;
import dr.evolution.coalescent.DemographicFunction;
import dr.evolution.coalescent.IntervalList;
import dr.evolution.coalescent.IntervalType;
import dr.evolution.tree.FlexibleNode;
import dr.evolution.tree.FlexibleTree;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evomodel.coalescent.DemographicModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.math.Binomial;
import dr.math.distributions.NormalGammaDistribution;
import dr.math.functionEval.GammaFunction;
import dr.xml.*;
import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.TDistributionImpl;

import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

/**
* General category model
*
* @author Matthew Hall
*/

public class GeneralCaseToCase extends CaseToCaseTreeLikelihood {

    public static final String GENERAL_CASE_TO_CASE = "generalCaseToCase";
    private Double[] timingLogLikelihoods;
    private Double[] storedTimingLogLikelihoods;


    public GeneralCaseToCase(PartitionedTreeModel virusTree, AbstractOutbreak caseData, String startingNetworkFileName,
                             Parameter infectionTimeBranchPositions, Parameter maxFirstInfToRoot)
            throws TaxonList.MissingTaxonException {
        this(virusTree, caseData, startingNetworkFileName, infectionTimeBranchPositions, null,
                maxFirstInfToRoot);
    }

    public GeneralCaseToCase(PartitionedTreeModel virusTree, AbstractOutbreak caseData, String startingNetworkFileName,
                             Parameter infectionTimeBranchPositions, Parameter infectiousTimePositions,
                             Parameter maxFirstInfToRoot)
            throws TaxonList.MissingTaxonException {
        super(GENERAL_CASE_TO_CASE, virusTree, caseData, infectionTimeBranchPositions, infectiousTimePositions,
                maxFirstInfToRoot);
        timingLogLikelihoods = new Double[noTips];
        storedTimingLogLikelihoods = new Double[noTips];

        prepareTree(startingNetworkFileName);

        prepareTimings();
    }

    public static double[] logOfAllValues(double[] values){
        double[] out = Arrays.copyOf(values, values.length);

        for(int i=0; i<values.length; i++){
            out[i] = Math.log(out[i]);
        }
        return out;
    }

    protected double calculateLogLikelihood(){

        // todo don't calculate this unless you have to

        double logL = 0;

        super.prepareTimings();

        int noInfectiousCategories = ((WithinCaseCategoryOutbreak) outbreak).getInfectiousCategoryCount();

        HashSet<String> infectiousCategories = ((WithinCaseCategoryOutbreak) outbreak).getInfectiousCategories();

        HashMap<String, NormalGammaDistribution> infectiousMap = ((WithinCaseCategoryOutbreak) outbreak).getInfectiousMap();

        HashMap<String, double[]> infectiousPriorParameters =
                new HashMap<String, double[]>();

        HashMap<String, Double> runningTotals = new HashMap<String, Double>();

        HashMap<String, Integer> counts = new HashMap<String, Integer>();

        for(String category: infectiousCategories){
            infectiousPriorParameters.put(category, Arrays.copyOf(infectiousMap.get(category).getParameters(), 4));
            runningTotals.put(category, 0.0);
            counts.put(category, 0);
        }

        boolean rootDone = false;
        HashSet<AbstractCase> remainingPossibleFirstInfections = new HashSet<AbstractCase>(outbreak.getCases());
        HashSet<AbstractCase> doneCases = new HashSet<AbstractCase>();




        // todo you really need to test whether order matters, if this thing ever works

        for(int i=0; i< outbreak.size(); i++){

            AbstractCase thisCase = outbreak.getCase(i);

            String category = ((WithinCaseCategoryOutbreak) outbreak).getInfectiousCategory(thisCase);

            double[] parameters = infectiousPriorParameters.get(category);

            Double runningTotal = runningTotals.get(category);
            Integer count = counts.get(category);

            double mu0 = infectiousMap.get(category).getParameters()[0];
            double lambda0 = infectiousMap.get(category).getParameters()[1];

            double oldMu = parameters[0];
            double oldLambda = parameters[1];
            double oldAlpha = parameters[2];
            double oldBeta = parameters[3];

            double infectiousPeriod = infectiousPeriods[i];

            NodeRef tip = treeModel.getNode(tipMap.get(thisCase));

            double startHeight = treeModel.getNodeHeight(tip);

            boolean maxFound = false;
            NodeRef currentNode = tip;
            double earliestInfectionHeight = startHeight;

            while(!maxFound){
                NodeRef parent = treeModel.getParent(currentNode);
                if(parent!=null){
                    AbstractCase parentCase = branchMap.get(parent.getNumber());
                    if(parentCase!=thisCase && doneCases.contains(parentCase)){
                        maxFound = true;
                        earliestInfectionHeight = treeModel.getNodeHeight(parent);
                    } else {
                        currentNode = parent;
                    }
                } else {
                    earliestInfectionHeight = treeModel.getNodeHeight(treeModel.getRoot()) +
                            maxFirstInfToRoot.getParameterValue(0);
                }
            }

            double maxInfPeriod = earliestInfectionHeight - startHeight;

            double minInfPeriod = 0;
            if(remainingPossibleFirstInfections.contains(thisCase) && remainingPossibleFirstInfections.size()==1){
                minInfPeriod = treeModel.getNodeHeight(treeModel.getRoot());
            }

            double numerator = tDistributionPDF(infectiousPeriod, oldMu, oldBeta*(oldLambda+1)/(oldAlpha*oldLambda),
                    2*oldAlpha);
            double denominator = 0;
            try{
                // watch this fall flat on its back...

                denominator = tDistributionCDF(maxInfPeriod, oldMu, oldBeta*(oldLambda+1)/(oldAlpha*oldLambda),
                        2*oldAlpha);
                if(minInfPeriod>0){
                    denominator -= tDistributionCDF(minInfPeriod, oldMu, oldBeta*(oldLambda+1)/(oldAlpha*oldLambda),
                            2*oldAlpha);
                }
            } catch (MathException e){
                throw new RuntimeException("Mathematical failure calculating the CDF of the T distribution");
            }

            logL += Math.log(numerator) - Math.log(denominator);

            runningTotal += infectiousPeriod;
            count += 1;

            double mean = runningTotal/count;


            double newMu = (mu0*lambda0 + runningTotal)/(lambda0 + count);
            double newLambda = oldLambda + 1;
            double newAlpha = oldAlpha + 0.5;
            double newBeta = oldBeta + oldLambda*Math.pow(mean-oldMu,2)/(2*(oldLambda+1));

            double[] newParameters = {newMu, newLambda, newAlpha, newBeta};
            infectiousPriorParameters.put(category, newParameters);
            doneCases.add(thisCase);
            remainingPossibleFirstInfections.remove(thisCase);
            for(AbstractCase descendant : getDescendants(thisCase)){
                remainingPossibleFirstInfections.remove(descendant);
            }



        }


        if(hasLatentPeriods){
            int noLatentCategories = ((WithinCaseCategoryOutbreak) outbreak).getLatentCategoryCount();

            HashSet<String> latentCategories = ((WithinCaseCategoryOutbreak) outbreak).getLatentCategories();

            HashMap<String, ArrayList<Double>> latentPeriodsByCategory = new HashMap<String, ArrayList<Double>>();

            for(String latentCategory : ((WithinCaseCategoryOutbreak) outbreak).getLatentCategories()){
                latentPeriodsByCategory.put(latentCategory, new ArrayList<Double>());
            }

            for(AbstractCase aCase : outbreak.getCases()){
                String category = ((WithinCaseCategoryOutbreak) outbreak).getLatentCategory(aCase);

                ArrayList<Double> correspondingList
                        = latentPeriodsByCategory.get(category);

                correspondingList.add(latentPeriods[outbreak.getCaseIndex(aCase)]);
            }

            for(String category: ((WithinCaseCategoryOutbreak) outbreak).getInfectiousCategories()){
                ArrayList<Double> latPeriodsInThisCategory = latentPeriodsByCategory.get(category);

                double count = (double)latPeriodsInThisCategory.size();

                NormalGammaDistribution prior = ((WithinCaseCategoryOutbreak) outbreak)
                        .getLatentCategoryPrior(category);

                double[] latPredictiveDistributionParameters=prior.getParameters();

                double mu_0 = latPredictiveDistributionParameters[0];
                double lambda_0 = latPredictiveDistributionParameters[1];
                double alpha_0 = latPredictiveDistributionParameters[2];
                double beta_0 = latPredictiveDistributionParameters[3];

                double lambda_n = lambda_0 + count;
                double alpha_n = alpha_0 + count/2;
                double sum = 0;
                for (Double latPeriod : latPeriodsInThisCategory) {
                    sum += latPeriod;
                }
                double mean = sum/count;

                double sumOfDifferences = 0;
                for (Double latPeriod : latPeriodsInThisCategory) {
                    sumOfDifferences += Math.pow(latPeriod-mean,2);
                }

                double mu_n = (lambda_0*mu_0 + sum)/(lambda_0 + count);
                double beta_n = beta_0 + 0.5*sumOfDifferences + lambda_0*count*Math.pow(mean-mu_0, 2)
                        /(2*(lambda_0+count));

                double priorPredictiveProbability
                        = GammaFunction.logGamma(alpha_n)
                        - GammaFunction.logGamma(alpha_0)
                        + alpha_0*Math.log(beta_0)
                        - alpha_n*Math.log(beta_n)
                        + 0.5*Math.log(lambda_0)
                        - 0.5*Math.log(lambda_n)
                        - (count/2)*Math.log(2*Math.PI);

                logL += priorPredictiveProbability;

                //todo log the parameters of the "posterior"
            }

        }


        likelihoodKnown = true;

        if(DEBUG){
            debugOutputTree("out.nex", true);
        }

        return logL;


    }

    private double tDistributionPDF(double x, double mu, double sigma, double df){
        TDistributionImpl distribution = new TDistributionImpl(df);
        return distribution.density((x-mu)/sigma);
    }

    private double tDistributionCDF(double x, double mu, double sigma, double df) throws MathException{
        TDistributionImpl distribution = new TDistributionImpl(df);
        return distribution.cumulativeProbability((x-mu)/sigma);
    }

    public void storeState(){
        super.storeState();
        storedTimingLogLikelihoods = Arrays.copyOf(timingLogLikelihoods, timingLogLikelihoods.length);
    }

    public void restoreState(){
        super.restoreState();
        timingLogLikelihoods = storedTimingLogLikelihoods;
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {

        // todo this is obviously wasteful at the moment

        super.handleModelChangedEvent(model, object, index);

        Arrays.fill(timingLogLikelihoods, null);
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {


        // todo this too is obviously wasteful at the moment

        super.handleVariableChangedEvent(variable, index, type);

        if(variable == infectionTimeBranchPositions || variable == infectiousTimePositions){

            Arrays.fill(timingLogLikelihoods, null);
        }
    }

    public void makeDirty(){
        super.makeDirty();
        Arrays.fill(timingLogLikelihoods, null);

    }

    // Tears the tree into small pieces. Indexes correspond to indexes in the outbreak.
    // todo Work out when components of this are unchanged after PT or TT moves



    public ArrayList<AbstractCase> postOrderTransmissionTreeTraversal(){
        return traverseTransmissionTree(branchMap.get(treeModel.getRoot().getNumber()));
    }

    private ArrayList<AbstractCase> traverseTransmissionTree(AbstractCase aCase){
        ArrayList<AbstractCase> out = new ArrayList<AbstractCase>();
        HashSet<AbstractCase> children = getInfectees(aCase);
        for(int i=0; i<getOutbreak().size(); i++){
            AbstractCase possibleChild = getOutbreak().getCase(i);
            // easiest way to maintain the set ordering of the outbreak?
            if(children.contains(possibleChild)){
                out.addAll(traverseTransmissionTree(possibleChild));
            }
        }
        out.add(aCase);
        return out;
    }

    private void copyPartitionToLittleTree(FlexibleTree littleTree, NodeRef oldNode, NodeRef newParent,
                                           AbstractCase partition){
        if(branchMap.get(oldNode.getNumber())==partition){
            if(treeModel.isExternal(oldNode)){
                NodeRef newTip = new FlexibleNode(new Taxon(treeModel.getNodeTaxon(oldNode).getId()));
                littleTree.addChild(newParent, newTip);
                littleTree.setBranchLength(newTip, treeModel.getBranchLength(oldNode));
            } else {
                NodeRef newChild = new FlexibleNode();
                littleTree.addChild(newParent, newChild);
                littleTree.setBranchLength(newChild, treeModel.getBranchLength(oldNode));
                for(int i=0; i<treeModel.getChildCount(oldNode); i++){
                    copyPartitionToLittleTree(littleTree, treeModel.getChild(oldNode, i), newChild, partition);
                }
            }
        } else {
            // we need a new tip
            NodeRef transmissionTip = new FlexibleNode(
                    new Taxon("Transmission_"+branchMap.get(oldNode.getNumber()).getName()));
            double parentTime = getNodeTime(treeModel.getParent(oldNode));
            double childTime = getInfectionTime(branchMap.get(oldNode.getNumber()));
            littleTree.addChild(newParent, transmissionTip);
            littleTree.setBranchLength(transmissionTip, childTime - parentTime);

        }
    }

    private class TreePlusRootBranchLength extends FlexibleTree {

        private double rootBranchLength;

        private TreePlusRootBranchLength(FlexibleTree tree, double rootBranchLength){
            super(tree);
            this.rootBranchLength = rootBranchLength;
        }

        private double getRootBranchLength(){
            return rootBranchLength;
        }

        private void setRootBranchLength(double rootBranchLength){
            this.rootBranchLength = rootBranchLength;
        }
    }

    private class MaxTMRCACoalescent extends Coalescent {

        private double maxHeight;

        private MaxTMRCACoalescent(Tree tree, DemographicModel demographicModel, double maxHeight){
            super(tree, demographicModel.getDemographicFunction());

            this.maxHeight = maxHeight;

        }

        public double calculateLogLikelihood() {
            return calculatePartitionTreeLogLikelihood(getIntervals(), getDemographicFunction(), 0, maxHeight);
        }

    }

    public static double calculatePartitionTreeLogLikelihood(IntervalList intervals,
                                                             DemographicFunction demographicFunction, double threshold,
                                                             double maxHeight) {

        double logL = 0.0;

        double startTime = -maxHeight;
        final int n = intervals.getIntervalCount();

        //TreeIntervals sets up a first zero-length interval with a lineage count of zero - skip this one

        for (int i = 0; i < n; i++) {

            // time zero corresponds to the date of first infection

            final double duration = intervals.getInterval(i);
            final double finishTime = startTime + duration;

            final double intervalArea = demographicFunction.getIntegral(startTime, finishTime);
            double normalisationArea = demographicFunction.getIntegral(startTime, 0);

            if( intervalArea == 0 && duration != 0 ) {
                return Double.NEGATIVE_INFINITY;
            }
            final int lineageCount = intervals.getLineageCount(i);

            if(lineageCount>=2){

                final double kChoose2 = Binomial.choose2(lineageCount);

                if (intervals.getIntervalType(i) == IntervalType.COALESCENT) {

                    logL += -kChoose2 * intervalArea;

                    final double demographicAtCoalPoint = demographicFunction.getDemographic(finishTime);

                    if( duration == 0.0 || demographicAtCoalPoint * (intervalArea/duration) >= threshold ) {
                        logL -= Math.log(demographicAtCoalPoint);
                    } else {
                        return Double.NEGATIVE_INFINITY;
                    }

                } else {
                    double numerator = Math.exp(-kChoose2 * intervalArea) - Math.exp(-kChoose2 * normalisationArea);
                    logL += Math.log(numerator);

                }

                // normalisation

                double logDenominator = Math.log1p(-Math.exp(-kChoose2 * normalisationArea));

                logL -= logDenominator;

            }

            startTime = finishTime;
        }

        return logL;
    }

    public void debugTreelet(Tree treelet, String fileName){
        try{
            FlexibleTree treeCopy = new FlexibleTree(treelet);
            for(int j=0; j<treeCopy.getNodeCount(); j++){
                FlexibleNode node = (FlexibleNode)treeCopy.getNode(j);
                node.setAttribute("Number", node.getNumber());
            }
            NexusExporter testTreesOut = new NexusExporter(new PrintStream(fileName));
            testTreesOut.exportTree(treeCopy);
        } catch (IOException ignored) {System.out.println("IOException");}
    }


    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {
        public static final String STARTING_NETWORK = "startingNetwork";
        public static final String INFECTION_TIMES = "infectionTimeBranchPositions";
        public static final String INFECTIOUS_TIMES = "infectiousTimePositions";
        public static final String MAX_FIRST_INF_TO_ROOT = "maxFirstInfToRoot";

        public String getParserName() {
            return GENERAL_CASE_TO_CASE;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            PartitionedTreeModel virusTree = (PartitionedTreeModel) xo.getChild(TreeModel.class);

            String startingNetworkFileName=null;

            if(xo.hasChildNamed(STARTING_NETWORK)){
                startingNetworkFileName = (String) xo.getElementFirstChild(STARTING_NETWORK);
            }

            AbstractOutbreak caseSet = (AbstractOutbreak) xo.getChild(AbstractOutbreak.class);

            CaseToCaseTreeLikelihood likelihood;

            Parameter infectionTimes = (Parameter) xo.getElementFirstChild(INFECTION_TIMES);

            Parameter infectiousTimes = xo.hasChildNamed(INFECTIOUS_TIMES)
                    ? (Parameter) xo.getElementFirstChild(INFECTIOUS_TIMES) : null;

            Parameter earliestFirstInfection = (Parameter) xo.getElementFirstChild(MAX_FIRST_INF_TO_ROOT);


            try {
                likelihood = new GeneralCaseToCase(virusTree, caseSet, startingNetworkFileName, infectionTimes,
                        infectiousTimes, earliestFirstInfection);
            } catch (TaxonList.MissingTaxonException e) {
                throw new XMLParseException(e.toString());
            }

            return likelihood;
        }

        public String getParserDescription() {
            return "This element provides a tree prior for a partitioned tree, with each partitioned tree generated" +
                    "by a coalescent process";
        }

        public Class getReturnType() {
            return GeneralCaseToCase.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                new ElementRule(PartitionedTreeModel.class, "The tree"),
                new ElementRule(WithinCaseCategoryOutbreak.class, "The set of outbreak"),
                new ElementRule("startingNetwork", String.class, "A CSV file containing a specified starting network",
                        true),
                new ElementRule(MAX_FIRST_INF_TO_ROOT, Parameter.class, "The maximum time from the first infection to" +
                        "the root node"),
                new ElementRule(INFECTION_TIMES, Parameter.class),
                new ElementRule(INFECTIOUS_TIMES, Parameter.class, "For each case, proportions of the time between " +
                        "infection and first event that requires infectiousness (further infection or cull)" +
                        "that has elapsed before infectiousness", true),
        };
    };
}
