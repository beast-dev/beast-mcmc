package dr.evomodel.epidemiology.casetocase;

import dr.evolution.util.TaxonList;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.xml.*;
import org.apache.commons.math.stat.descriptive.DescriptiveStatistics;

/**
 * Created by mhall on 17/12/2013.
 */
public class JeffreysCaseToCase extends CaseToCaseTreeLikelihood  {

    public static final String JEFFREYS_CASE_TO_CASE = "jeffreysCaseToCase";

    public JeffreysCaseToCase(PartitionedTreeModel virusTree, AbstractOutbreak caseData, String startingNetworkFileName,
                                    Parameter infectionTimeBranchPositions, Parameter infectiousTimePositions,
                                    Parameter maxFirstInfToRoot) throws TaxonList.MissingTaxonException{
        super(virusTree, caseData, infectionTimeBranchPositions, infectiousTimePositions,
                maxFirstInfToRoot);

        prepareTree(startingNetworkFileName);
    }


    protected double calculateLogLikelihood() {

        // todo implement categories

        super.prepareTimings();

        if(!hasLatentPeriods){
            return isAllowed() ? -Math.log(getLogInfectiousPeriodSD()) : Double.NEGATIVE_INFINITY;
        } else {
            double infL = -Math.log(getLogInfectiousPeriodSD());
            double latL = -Math.log(getLogLatentPeriodSD());
            return isAllowed() ? infL + latL : Double.NEGATIVE_INFINITY;
        }
    }


    public double getInfectiousPeriodSD(){
        DescriptiveStatistics stats = new DescriptiveStatistics(convertArray(getInfectiousPeriods(true)));
        return stats.getStandardDeviation();
    }

    public double getLogInfectiousPeriodSD(){
        if(infectiousPeriods==null){
            infectiousPeriods = getInfectiousPeriods(true);
        }
        DescriptiveStatistics stats = new DescriptiveStatistics();
        for (double infectiousPeriod : infectiousPeriods) {
            stats.addValue(Math.log(infectiousPeriod));
        }
        return stats.getStandardDeviation();
    }

    public double getLatentPeriodSD(){
        if(latentPeriods==null){
            latentPeriods = getLatentPeriods(true);
        }
        DescriptiveStatistics stats = new DescriptiveStatistics(convertArray(getLatentPeriods(true)));
        return stats.getStandardDeviation();
    }

    public double getLogLatentPeriodSD(){
        if(latentPeriods==null){
            latentPeriods = getLatentPeriods(true);
        }
        DescriptiveStatistics stats = new DescriptiveStatistics();
        for (double latentPeriod : latentPeriods) {
            stats.addValue(Math.log(latentPeriod));
        }
        return stats.getStandardDeviation();
    }

    // todo This is already pretty fast, but changeMap could be overriden.

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {
        public static final String STARTING_NETWORK = "startingNetwork";
        public static final String INFECTION_TIMES = "infectionTimeBranchPositions";
        public static final String INFECTIOUS_TIMES = "infectiousTimePositions";
        public static final String MAX_FIRST_INF_TO_ROOT = "maxFirstInfToRoot";

        public String getParserName() {
            return JEFFREYS_CASE_TO_CASE;
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
                likelihood = new JeffreysCaseToCase(virusTree, caseSet, startingNetworkFileName, infectionTimes,
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
            return WithinCaseCoalescent.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                new ElementRule(PartitionedTreeModel.class, "The tree"),
                new ElementRule(JeffreysCategoryOutbreak.class, "The set of cases"),
                new ElementRule("startingNetwork", String.class, "A CSV file containing a specified starting network",
                        true),
                new ElementRule(MAX_FIRST_INF_TO_ROOT, Parameter.class, "The maximum time from the first infection to" +
                        "the root node"),
                new ElementRule(INFECTION_TIMES, Parameter.class),
                new ElementRule(INFECTIOUS_TIMES, Parameter.class, "For each case, proportions of the time between " +
                        "infection and first event that requires infectiousness (further infection or cull)" +
                        "that has elapsed before infectiousness", true)
        };
    };



}
