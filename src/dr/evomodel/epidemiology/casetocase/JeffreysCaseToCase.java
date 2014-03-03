package dr.evomodel.epidemiology.casetocase;

import dr.evolution.util.TaxonList;
import dr.evomodel.epidemiology.casetocase.periodpriors.AbstractPeriodPriorDistribution;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.ArrayList;
import java.util.HashMap;

/**
 * Created by mhall on 17/12/2013.
 *
 * todo there's a better OO structure out there - WCC could be an extension of this class
 */
public class JeffreysCaseToCase extends CaseToCaseTreeLikelihood  {

    private double infectiousPeriodsLogLikelihood;
    private double storedInfectiousPeriodsLogLikelihood;
    private double latentPeriodsLogLikelihood;
    private double storedLatentPeriodsLogLikelihood;

    public static final String JEFFREYS_CASE_TO_CASE = "jeffreysCaseToCase";

    public JeffreysCaseToCase(PartitionedTreeModel virusTree, AbstractOutbreak caseData, String startingNetworkFileName,
                                    Parameter infectionTimeBranchPositions, Parameter infectiousTimePositions,
                                    Parameter maxFirstInfToRoot)
            throws TaxonList.MissingTaxonException{
        super(virusTree, caseData, infectionTimeBranchPositions, infectiousTimePositions,
                maxFirstInfToRoot);

        prepareTree(startingNetworkFileName);
    }


    protected double calculateLogLikelihood() {

        // you shouldn't need to do this, because C2CTransL will already have done it

        // super.prepareTimings();

        HashMap<String, ArrayList<Double>> infectiousPeriodsByCategory
                = new HashMap<String, ArrayList<Double>>();

        // todo do this only once? Using indexes?

        for(AbstractCase aCase : outbreak.getCases()){

            String category = ((CategoryOutbreak) outbreak).getInfectiousCategory(aCase);

            if(!infectiousPeriodsByCategory.keySet().contains(category)){
                infectiousPeriodsByCategory.put(category, new ArrayList<Double>());
            }

            ArrayList<Double> correspondingList
                    = infectiousPeriodsByCategory.get(category);

            correspondingList.add(getInfectiousPeriod(aCase));
        }

        infectiousPeriodsLogLikelihood = 0;

        for(String category : ((CategoryOutbreak) outbreak).getInfectiousCategories()){

            // todo this inevitably should be an improper distribution, and a warning if it isn't is in order

            Double[] infPeriodsInThisCategory = infectiousPeriodsByCategory.get(category)
                    .toArray(new Double[infectiousPeriodsByCategory.size()]);

            AbstractPeriodPriorDistribution hyperprior = ((CategoryOutbreak) outbreak)
                    .getInfectiousCategoryPrior(category);

            double[] values = new double[infPeriodsInThisCategory.length];

            for(int i=0; i<infPeriodsInThisCategory.length; i++){
                values[i] = infPeriodsInThisCategory[i];
            }

            infectiousPeriodsLogLikelihood += hyperprior.getLogLikelihood(values);

        }


        if(hasLatentPeriods){

            HashMap<String, ArrayList<Double>> latentPeriodsByCategory
                    = new HashMap<String, ArrayList<Double>>();

            // todo do this only once?

            for(AbstractCase aCase : outbreak.getCases()){

                String category = ((CategoryOutbreak) outbreak).getLatentCategory(aCase);

                if(!latentPeriodsByCategory.keySet().contains(category)){
                    latentPeriodsByCategory.put(category, new ArrayList<Double>());
                }

                ArrayList<Double> correspondingList
                        = latentPeriodsByCategory.get(category);

                correspondingList.add(getLatentPeriod(aCase));
            }


            latentPeriodsLogLikelihood = 0;

            for(String category : ((CategoryOutbreak) outbreak).getLatentCategories()){

                Double[] latPeriodsInThisCategory = latentPeriodsByCategory.get(category)
                        .toArray(new Double[latentPeriodsByCategory.size()]);

                AbstractPeriodPriorDistribution hyperprior = ((CategoryOutbreak) outbreak)
                        .getLatentCategoryPrior(category);

                double[] values = new double[latPeriodsInThisCategory.length];

                for(int i=0; i<latPeriodsInThisCategory.length; i++){
                    values[i] = latPeriodsInThisCategory[i];
                }

                latentPeriodsLogLikelihood += hyperprior.getLogLikelihood(values);

            }

        }

        return infectiousPeriodsLogLikelihood + latentPeriodsLogLikelihood;

    }

    public void storeState(){
        super.storeState();
        storedInfectiousPeriodsLogLikelihood = infectiousPeriodsLogLikelihood;
        if(hasLatentPeriods){
            storedLatentPeriodsLogLikelihood = latentPeriodsLogLikelihood;
        }
    }

    public void restoreState(){
        super.restoreState();
        infectiousPeriodsLogLikelihood = storedInfectiousPeriodsLogLikelihood;
        if(hasLatentPeriods){
            latentPeriodsLogLikelihood = storedLatentPeriodsLogLikelihood;
        }
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
            return JEFFREYS_CASE_TO_CASE;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            PartitionedTreeModel virusTree = (PartitionedTreeModel) xo.getChild(TreeModel.class);

            String startingNetworkFileName = null;

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
            return "This element provides the likelihood of a partitioned tree.";
        }

        public Class getReturnType() {
            return WithinCaseCoalescent.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                new ElementRule(PartitionedTreeModel.class, "The tree"),
                new ElementRule(CategoryOutbreak.class, "The set of outbreak"),
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




