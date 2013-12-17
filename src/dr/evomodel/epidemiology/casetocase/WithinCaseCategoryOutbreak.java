package dr.evomodel.epidemiology.casetocase;

import dr.evolution.util.Date;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.inference.distribution.ParametricDistributionModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.ProductStatistic;
import dr.inference.model.Variable;
import dr.math.IntegrableUnivariateFunction;
import dr.math.Integral;
import dr.math.RiemannApproximation;
import dr.math.UnivariateFunction;
import dr.xml.*;

import java.util.ArrayList;

/**
 * Outbreak class for within-case coalescent.
 *
 * Each case belongs to an infectious (and latent) category which corresponds to one of a list of probability
 * distributions (most likely gamma or exponential) for the length of the infectious (latent) period. The XML rules for
 * the outbreak class ask for at least one ParametricDistributionModel.
 * Assignment of cases to distributions should be handled in whatever script or GUI writes the XML.
 *
 * Intended for situations where no data on infection times exists.
 *
 * User: Matthew Hall
 * Date: 17/12/2013
 */

public class WithinCaseCategoryOutbreak extends AbstractOutbreak {

    public static final String WITHIN_CASE_CATEGORY_OUTBREAK = "withinCaseCategoryOutbreak";
    private double[][] distances;

    public WithinCaseCategoryOutbreak(String name, Taxa taxa, boolean hasGeography, boolean hasLatentPeriods){
        super(name, taxa, hasLatentPeriods, hasGeography);
        cases = new ArrayList<AbstractCase>();
    }


    public WithinCaseCategoryOutbreak(String name, Taxa taxa, boolean hasGeography, boolean hasLatentPeriods,
                                      ArrayList<AbstractCase> cases){
        this(name, taxa, hasGeography, hasLatentPeriods);
        this.cases.addAll(cases);
    }


    private void addCase(String caseID, Date examDate, Date cullDate, ParametricDistributionModel infectiousDist,
                         ParametricDistributionModel latentDist, Parameter coords, Taxa associatedTaxa){
        WithinCaseCategoryCase thisCase;

        if(latentDist==null){
            thisCase =  new WithinCaseCategoryCase(caseID, examDate, cullDate, infectiousDist, coords, associatedTaxa);
        } else {
            thisCase =
                    new WithinCaseCategoryCase(caseID, examDate, cullDate, infectiousDist, latentDist, coords, associatedTaxa);
        }
        cases.add(thisCase);
        addModel(thisCase);
    }



    public double getDistance(AbstractCase a, AbstractCase b) {
        if(distances==null){
            throw new RuntimeException("Distance matrix has not been initialised");
        }
        return distances[getCaseIndex(a)][getCaseIndex(b)];
    }

    private void buildDistanceMatrix(){
        distances = new double[cases.size()][cases.size()];

        for(int i=0; i<cases.size(); i++){
            for(int j=0; j<cases.size(); j++){
                distances[i][j]=SpatialKernel.EuclideanDistance(getCase(i).getCoords(),getCase(j).getCoords());
            }
        }
    }

    public double infectedAtGivenBefore(int caseIndex, double infected, double before){
        WithinCaseCategoryCase thisCase = (WithinCaseCategoryCase)cases.get(caseIndex);
        return thisCase.infectedAtGivenBefore(infected, before);
    }

    public double infectedAtGivenBeforeInfectiousAtGivenBefore(int caseIndex, double infected,
                                                               double infectedBefore, double infectious,
                                                               double infectiousBefore){
        WithinCaseCategoryCase thisCase = (WithinCaseCategoryCase)cases.get(caseIndex);
        return thisCase.infectedAtGivenBeforeInfectiousAtGivenBefore(infected, infectedBefore, infectious,
                infectiousBefore);
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        fireModelChanged();
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        fireModelChanged();
    }

    protected void storeState() {
        //nothing to do
    }

    protected void restoreState() {
        //nothing to do
    }

    protected void acceptState() {
        //nothing to do
    }

    private class WithinCaseCategoryCase extends AbstractCase{

        public static final String WITHIN_CASE_CATEGORY_CASE = "withinCaseCategoryCase";
        private ParametricDistributionModel infectiousPeriodDistribution;
        private ParametricDistributionModel latentPeriodDistribution;
        private Parameter coords;

        private WithinCaseCategoryCase(String name, String caseID, Date examDate, Date cullDate,
                                       ParametricDistributionModel infectiousDist, Parameter coords,
                                       Taxa associatedTaxa){
            super(name);
            this.caseID = caseID;
            this.examDate = examDate;
            endOfInfectiousDate = cullDate;
            this.associatedTaxa = associatedTaxa;
            this.coords = coords;
            infectiousPeriodDistribution = infectiousDist;
            this.addModel(infectiousPeriodDistribution);
            latentPeriodDistribution = null;
        }


        private WithinCaseCategoryCase(String name, String caseID, Date examDate, Date cullDate,
                                       ParametricDistributionModel infectiousDist,
                                       ParametricDistributionModel latentDist, Parameter coords,
                                       Taxa associatedTaxa){
            this(name, caseID, examDate, cullDate, infectiousDist, coords, associatedTaxa);
            this.latentPeriodDistribution = latentDist;
            this.addModel(latentPeriodDistribution);
        }


        private WithinCaseCategoryCase(String caseID, Date examDate, Date cullDate,
                                       ParametricDistributionModel infectiousDist, Parameter coords,
                                       Taxa associatedTaxa){
            this(WITHIN_CASE_CATEGORY_CASE, caseID, examDate, cullDate, infectiousDist, coords, associatedTaxa);
        }


        private WithinCaseCategoryCase(String caseID, Date examDate, Date cullDate,
                                       ParametricDistributionModel infectiousDist,
                                       ParametricDistributionModel latentDist, Parameter coords, Taxa associatedTaxa){
            this(WITHIN_CASE_CATEGORY_CASE, caseID, examDate, cullDate, infectiousDist, latentDist, coords,
                    associatedTaxa);
        }

        public double infectedAt(double infected){
            if(culledYet(infected)){
                return 0;
            } else {
                return infectiousPeriodDistribution.pdf(endOfInfectiousDate.getTimeValue()-infected);
            }
        }

        public boolean culledYet(double time) {
            return time>endOfInfectiousDate.getTimeValue();
        }

        public boolean examinedYet(double time) {
            return time>examDate.getTimeValue();
        }


        //probability that this case was infected at "infected" given that it was infected before "after" and
        // there are no latent periods

        public double infectedAtGivenBefore(double infected, double before){
            if(infected>before){
                return 0;
            } else {
                if(latentPeriodDistribution == null){
                    return infectedAt(infected)
                            /(1-infectiousPeriodDistribution.cdf(endOfInfectiousDate.getTimeValue()-before));
                } else {
                    throw new RuntimeException("Calling the probability function for a model without latent" +
                            "periods, but there is one");
                }
            }
        }

        //probability that this case was infected at "infected" given that it was infected after "infectedAfter" and
        // became infectious at "infectious" given that it was infectious after "infectiousAfter"

        public double infectedAtGivenBeforeInfectiousAtGivenBefore(double infected, double infectedBefore,
                                                                   double infectious, double infectiousBefore){
            if(infected>infectedBefore || infectious>infectiousBefore || infectious<infected){
                return 0;
            } else {
                if(latentPeriodDistribution != null){
                    // probability of infection date given known infectious date & cull date of known parent
                    double eoi = getCullTime();
                    double latentPeriod = infectious-infected;
                    double minLatentPeriod = Math.min(infectious-infectedBefore, 0);
                    double out = latentPeriodDistribution.pdf(latentPeriod)
                            /(1-latentPeriodDistribution.cdf(minLatentPeriod));
                    // probability of infectious date given known cull date & known infection dates of known
                    // children
                    double infectiousPeriod = eoi - infectious;
                    double minInfectiousPeriod = Math.min(eoi - infectiousBefore, 0);
                    out *= infectiousPeriodDistribution.pdf(infectiousPeriod)
                            /(1-infectiousPeriodDistribution.cdf(minInfectiousPeriod));
                    return out;

                } else {
                    throw new RuntimeException("Calling the probability function for a model with latent" +
                            "periods, but we have none");
                }
            }
        }

        protected void handleModelChangedEvent(Model model, Object object, int index) {
            // @todo to have all the cases listening seems excessive and I'm not sure it's necessary - maybe only the outbreak need listen
            fireModelChanged();
        }

        protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
            fireModelChanged();
        }

        protected void storeState() {
            // nothing to do?
        }

        protected void restoreState() {
            // nothing to do?
        }

        protected void acceptState() {
            //nothing to do?
        }

        public double[] getCoords() {
            return new double[]{coords.getParameterValue(0), coords.getParameterValue(1)};
        }
    }

// for integrating out infectiousness dates, if this is preferred:

    private class CombinedPeriodFunction implements IntegrableUnivariateFunction {

        private ParametricDistributionModel infectious;
        private ParametricDistributionModel latent;
        private Integral numIntergrator;
        private PdfByPdf pdfByPdf;
        private CdfByPdf cdfByPdf;

        private CombinedPeriodFunction(ParametricDistributionModel infectious, ParametricDistributionModel latent,
                                       int numSteps){
            this.infectious = infectious;
            this.latent = latent;
            this.numIntergrator = new RiemannApproximation(numSteps);
            cdfByPdf = new CdfByPdf(1);
            pdfByPdf = new PdfByPdf(1);
        }

        public double evaluateIntegral(double a, double b) {
            cdfByPdf.setTotal(b);
            double out = numIntergrator.integrate(cdfByPdf, 0, b);
            if(a>0){
                cdfByPdf.setTotal(a);
                out -= numIntergrator.integrate(cdfByPdf, 0, a);
            }
            return out;
        }

        public double evaluateIntegral(double a, double b, double maxLatent){
            double out;
            cdfByPdf.setTotal(b);
            if(maxLatent > b){
                out = numIntergrator.integrate(cdfByPdf, 0, b);
                out /= latent.cdf(maxLatent);
            } else {
                out = numIntergrator.integrate(cdfByPdf, 0, maxLatent);
                out /= latent.cdf(maxLatent);
            }
            if(a>0){
                cdfByPdf.setTotal(a);
                if(maxLatent > a){
                    out -= numIntergrator.integrate(cdfByPdf, 0, a)/latent.cdf(maxLatent);
                } else {
                    out -= numIntergrator.integrate(cdfByPdf, 0, maxLatent)/latent.cdf(maxLatent);
                }

            }
            return out;
        }

        public double evaluate(double argument) {
            pdfByPdf.setTotal(argument);
            return numIntergrator.integrate(pdfByPdf, 0, argument);
        }

        public double evaluate(double argument, double maxLatent){
            if(maxLatent>argument){
                return evaluate(argument)/latent.cdf(maxLatent);
            }
            pdfByPdf.setTotal(argument);
            return numIntergrator.integrate(pdfByPdf, 0, maxLatent)/latent.cdf(maxLatent);
        }


        public double getLowerBound() {
            return 0;
        }

        public double getUpperBound() {
            return Double.POSITIVE_INFINITY;
        }

        private class PdfByPdf implements UnivariateFunction {

            double total;

            private PdfByPdf(double total){
                this.total = total;
            }

            public double evaluate(double argument) {
                return infectious.pdf(total-argument)*latent.pdf(argument);
            }

            public double getLowerBound() {
                return 0;
            }

            public double getUpperBound() {
                return Double.POSITIVE_INFINITY;
            }

            public void setTotal(double total){
                this.total = total;
            }

            public double getTotal(){
                return total;
            }

        }

        private class CdfByPdf implements UnivariateFunction{

            double total;

            private CdfByPdf(double total){
                this.total = total;
            }

            public double evaluate(double argument) {
                return infectious.cdf(total-argument)*latent.pdf(argument);
            }

            public double getLowerBound() {
                return 0;
            }

            public double getUpperBound() {
                return Double.POSITIVE_INFINITY;
            }

            public void setTotal(double total){
                this.total = total;
            }

            public double getTotal(){
                return total;
            }

        }
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        //for the outbreak

        public static final String HAS_GEOGRAPHY = "hasGeography";
        public static final String HAS_LATENT_PERIODS = "hasLatentPeriods";
        public static final String INFECTIOUS_PERIOD_DISTRIBUTIONS = "infectiousPeriodDistributions";
        public static final String LATENT_PERIOD_DISTRIBUTIONS = "latentPeriodDistributions";

        //for the cases

        public static final String CASE_ID = "caseID";
        public static final String CULL_DAY = "cullDay";
        public static final String EXAMINATION_DAY = "examinationDay";
        public static final String COORDINATES = "spatialCoordinates";
        public static final String INFECTION_TIME_BRANCH_POSITION = "infectionTimeBranchPosition";
        public static final String INFECTIOUS_TIME_POSITION = "infectiousTimePosition";
        public static final String INFECTIOUS_PERIOD_DISTRIBUTION = "infectiousPeriodDistribution";
        public static final String LATENT_PERIOD_DISTRIBUTION = "latentPeriodDistribution";

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            final boolean hasGeography = xo.hasAttribute(HAS_GEOGRAPHY) && xo.getBooleanAttribute(HAS_GEOGRAPHY);
            final boolean hasLatentPeriods = xo.hasChildNamed(LATENT_PERIOD_DISTRIBUTIONS);
            final Taxa taxa = (Taxa) xo.getChild(Taxa.class);
            WithinCaseCategoryOutbreak cases = new WithinCaseCategoryOutbreak(null, taxa, hasGeography,
                    hasLatentPeriods);
            for(int i=0; i<xo.getChildCount(); i++){
                Object cxo = xo.getChild(i);
                if(cxo instanceof XMLObject && ((XMLObject)cxo).getName()
                        .equals(WithinCaseCategoryCase.WITHIN_CASE_CATEGORY_CASE)){
                    parseCase((XMLObject)cxo, cases, hasLatentPeriods);
                }
            }
            cases.buildDistanceMatrix();
            return cases;
        }

        public void parseCase(XMLObject xo, WithinCaseCategoryOutbreak outbreak, boolean expectLatentPeriods)
                throws XMLParseException {
            String farmID = (String) xo.getAttribute(CASE_ID);
            final Date cullDate = (Date) xo.getElementFirstChild(CULL_DAY);
            final Date examDate = (Date) xo.getElementFirstChild(EXAMINATION_DAY);
            final ParametricDistributionModel infectiousDist =
                    (ParametricDistributionModel)xo.getElementFirstChild(INFECTIOUS_PERIOD_DISTRIBUTION);
            ParametricDistributionModel latentDist = null;
            if(xo.hasChildNamed(LATENT_PERIOD_DISTRIBUTION)){
                latentDist = (ParametricDistributionModel)xo.getElementFirstChild(LATENT_PERIOD_DISTRIBUTION);
            } else if(expectLatentPeriods){
                throw new XMLParseException("Case "+farmID+" not assigned a latent periods distribution");
            }
            if(expectLatentPeriods && !xo.hasChildNamed(INFECTIOUS_TIME_POSITION)){
                throw new XMLParseException("Latent periods specified, but case "+farmID+" not assigned a time of " +
                        "infectiousness");
            }

            final Parameter coords = xo.hasChildNamed(COORDINATES) ?
                    (Parameter) xo.getElementFirstChild(COORDINATES) : null;
            Taxa taxa = new Taxa();
            for(int i=0; i<xo.getChildCount(); i++){
                if(xo.getChild(i) instanceof Taxon){
                    taxa.addTaxon((Taxon)xo.getChild(i));
                }
            }
            outbreak.addCase(farmID, examDate, cullDate, infectiousDist, latentDist, coords, taxa);
        }



        public String getParserDescription(){
            return "Parses a set of 'category' farm cases and the information that they all share";
        }

        public Class getReturnType(){
            return SimpleOutbreak.class;
        }

        public String getParserName(){
            return WITHIN_CASE_CATEGORY_OUTBREAK;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] caseRules = {
                new StringAttributeRule(CASE_ID, "The unique identifier for this farm"),
                new ElementRule(CULL_DAY, Date.class, "The date this farm was culled", false),
                new ElementRule(EXAMINATION_DAY, Date.class, "The date this farm was examined", false),
                new ElementRule(Taxon.class, 0, Integer.MAX_VALUE),
                new ElementRule(INFECTIOUS_PERIOD_DISTRIBUTION, ParametricDistributionModel.class,
                        "The probability distribution from which the infectious period of this case is drawn"),
                new ElementRule(LATENT_PERIOD_DISTRIBUTION, ParametricDistributionModel.class,
                        "The probability distribution from which the latent period of this case is drawn", true),
                new ElementRule(INFECTION_TIME_BRANCH_POSITION, Parameter.class, "The exact position on the branch" +
                        " along which the infection of this case occurs that it actually does occur"),
                new ElementRule(INFECTIOUS_TIME_POSITION, Parameter.class, "Parameter taking a value between 0 and" +
                        "1, indicating when from infection (0) to first caused infection (or cull if the cases" +
                        "causes no infections) (1) the case became infectious", true),
                new ElementRule(COORDINATES, Parameter.class, "The spatial coordinates of this case", true)
        };

        private final XMLSyntaxRule[] rules = {
                new ElementRule(ProductStatistic.class, 0,2),
                new ElementRule(WithinCaseCategoryCase.WITHIN_CASE_CATEGORY_CASE, caseRules, 1, Integer.MAX_VALUE),
                new ElementRule(Taxa.class),
                new ElementRule(INFECTIOUS_PERIOD_DISTRIBUTIONS, ParametricDistributionModel.class,
                        "One or more probability distributions for the infectious periods of cases in the oubreak", 1,
                        Integer.MAX_VALUE),
                new ElementRule(LATENT_PERIOD_DISTRIBUTIONS, ParametricDistributionModel.class,
                        "Zero or more probability distributions for the latent periods of cases in the oubreak", 0,
                        Integer.MAX_VALUE),
                AttributeRule.newBooleanRule(HAS_GEOGRAPHY, true)
        };
    };


}
