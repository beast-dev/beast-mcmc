package dr.evomodel.epidemiology.casetocase;

import dr.evolution.util.Date;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evolution.util.Units;
import dr.inference.distribution.ParametricDistributionModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.ProductStatistic;
import dr.inference.model.Variable;
import dr.math.RiemannApproximation;
import dr.xml.*;

import java.util.ArrayList;

/**
 * Each case will belong to an infectious category (and in future could have a latent category as well) which
 * corresponds to one of a list of probability distributions (most likely gamma or exponential) for the length of the
 * infectious period. The XML rules for the outbreak class ask for at least one ParametricDistributionModel.
 * Assignment of cases to distributions should be handled in whatever script or GUI writes the XML.
 *
 * Intended for situations where no data on infection times exists.
 *
 * User: Matthew Hall
 * Date: 20/08/2013
 */

public class CategoryOutbreak extends AbstractOutbreak {

    public static final String CATEGORY_OUTBREAK = "categoryOutbreak";
    private RiemannApproximation integrator;

    public CategoryOutbreak(String name, Taxa taxa, boolean hasGeography, Parameter riemannSampleSize){
        super(name, taxa);
        integrator = new RiemannApproximation((int)riemannSampleSize.getParameterValue(0));
        cases = new ArrayList<AbstractCase>();
        hasLatentPeriods = false;
        this.hasGeography = hasGeography;
    }

    public CategoryOutbreak(String name, Taxa taxa, boolean hasGeography, Parameter riemannSampleSize,
                            ArrayList<AbstractCase> cases){
        this(name, taxa, hasGeography, riemannSampleSize);
        this.cases.addAll(cases);
    }

    public CategoryOutbreak(Taxa taxa, boolean hasGeography, Parameter riemannSampleSize){
        this(CATEGORY_OUTBREAK, taxa, hasGeography, riemannSampleSize);
    }

    public CategoryOutbreak(Taxa taxa, boolean hasGeography, Parameter riemannSampleSize,
                            ArrayList<AbstractCase> cases){
        this(CATEGORY_OUTBREAK, taxa, hasGeography, riemannSampleSize, cases);
    }

    private void addCase(String caseID, Date examDate, Date cullDate, ParametricDistributionModel infectiousDist,
                         Parameter coords, Taxa associatedTaxa){
        CategoryCase thisCase = new CategoryCase(caseID, examDate, cullDate, infectiousDist, coords, associatedTaxa);
        cases.add(thisCase);
        addModel(thisCase);
    }

    public double getDistance(AbstractCase a, AbstractCase b) {
        return SpatialKernel.EuclideanDistance(a.getCoords(), b.getCoords());
    }

    // in all of the following infectiousness of the parent is assumed because there is no latent period, so Y is only
    // used to determine whether it was culled

    @Override
    public double probXInfectedByYAtTimeT(AbstractCase X, AbstractCase Y, double T) {
        if(Y.culledYet(T)){
            return 0;
        } else {
            return probXInfectedAtTimeT(X, T);
        }
    }

    @Override
    public double logProbXInfectedByYAtTimeT(AbstractCase X, AbstractCase Y, double T) {
        if(Y.culledYet(T)){
            return Double.NEGATIVE_INFINITY;
        } else {
            return logProbXInfectedAtTimeT(X, T);
        }
    }

    @Override
    public double probXInfectedByYBetweenTandU(AbstractCase X, AbstractCase Y, double T, double U) {
        if(Y.culledYet(T)){
            return 0;
        } else {
            double latestInfectionDate = Math.min(U, Y.getCullDate().getTimeValue());
            return probXInfectedBetweenTandU(X, T, latestInfectionDate);
        }
    }

    @Override
    public double logProbXInfectedByYBetweenTandU(AbstractCase X, AbstractCase Y, double T, double U) {
        if(Y.culledYet(T)){
            return Double.NEGATIVE_INFINITY;
        } else {
            double latestInfectionDate = Math.min(U, Y.getCullDate().getTimeValue());
            return logProbXInfectedBetweenTandU(X, T, latestInfectionDate);
        }
    }

    @Override
    public double probXInfectiousByTimeT(AbstractCase X, double T) {
        return ((CategoryCase)X).infectedBy(T);
    }

    @Override
    public double logProbXInfectiousByTimeT(AbstractCase X, double T) {
        return Math.log(probXInfectiousByTimeT(X, T));
    }

    @Override
    public double probXInfectedAtTimeT(AbstractCase X, double T) {
        return ((CategoryCase)X).infectedAt(T);
    }

    @Override
    public double logProbXInfectedAtTimeT(AbstractCase X, double T) {
        return Math.log(probXInfectedAtTimeT(X,T));
    }

    @Override
    public double probXInfectedBetweenTandU(AbstractCase X, double T, double U) {
        return ((CategoryCase)X).infectedBetween(T,U);
    }

    @Override
    public double logProbXInfectedBetweenTandU(AbstractCase X, double T, double U) {
        return Math.log(probXInfectedBetweenTandU(X, T, U));
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
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

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        fireModelChanged();
    }

    private class CategoryCase extends AbstractCase{

        public static final String CATEGORY_CASE = "categoryCase";
        private final Parameter coords;
        private ParametricDistributionModel infectiousPeriodDistribution;

        private CategoryCase(String name, String caseID, Date examDate, Date cullDate,
                             ParametricDistributionModel infectiousDist, Parameter coords, Taxa associatedTaxa){
            super(name);
            this.caseID = caseID;
            this.examDate = examDate;
            endOfInfectiousDate = cullDate;
            this.associatedTaxa = associatedTaxa;
            this.coords = coords;
            infectiousPeriodDistribution = infectiousDist;
            this.addModel(infectiousPeriodDistribution);
        }

        private CategoryCase(String caseID, Date examDate, Date cullDate, ParametricDistributionModel infectiousDist,
                             Parameter coords, Taxa associatedTaxa){
            this(CATEGORY_CASE, caseID, examDate, cullDate, infectiousDist, coords, associatedTaxa);
        }

        public Date getLatestPossibleInfectionDate() {
            Double doubleDate = examDate.getTimeValue();
            return Date.createTimeSinceOrigin(doubleDate, Units.Type.DAYS, examDate.getOrigin());
        }

        public double infectedAt(double infected){
            if(culledYet(infected)){
                return 0;
            } else {
                return infectiousPeriodDistribution.pdf(endOfInfectiousDate.getTimeValue()-infected);
            }
        }

        public double infectedBetween(double start, double end){
            if(culledYet(start)){
                return 0;
            } else {
                double endPoint = end<endOfInfectiousDate.getTimeValue() ? end : endOfInfectiousDate.getTimeValue();
                return infectiousPeriodDistribution.cdf(endOfInfectiousDate.getTimeValue()-start)
                        - infectiousPeriodDistribution.cdf(endOfInfectiousDate.getTimeValue()-endPoint);
            }

        }

        public double infectedBy(double time){
            if(culledYet(time)){
                return 1;
            } else {
                return 1 - infectiousPeriodDistribution.cdf(endOfInfectiousDate.getTimeValue()-time);
            }
        }

        public double infectedAfter(double time){
            return 1 - infectedBy(time);
        }

        public boolean culledYet(double time) {
            return time>endOfInfectiousDate.getTimeValue();
        }

        public boolean examinedYet(double time) {
            return time>examDate.getTimeValue();
        }

        public double[] getCoords() {
            return new double[]{coords.getParameterValue(0), coords.getParameterValue(1)};
        }

        protected void handleModelChangedEvent(Model model, Object object, int index) {
            fireModelChanged();
        }

        protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
            fireModelChanged();
        }

        protected void storeState() {
            // nothing to do
        }

        protected void restoreState() {
            // nothing to do
        }

        protected void acceptState() {
            //nothing to do
        }

        public Date getExamDate() {
            return examDate;
        }

        public ParametricDistributionModel getInfectiousPeriodDistribution(){
            return infectiousPeriodDistribution;
        }
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        //for the outbreak

        public static final String RIEMANN_SAMPLE_SIZE = "riemannSampleSize";
        public static final String HAS_GEOGRAPHY = "hasGeography";
        public static final String INFECTIOUS_PERIOD_DISTRIBUTIONS = "infectiousPeriodDistributions";

        //for the cases

        public static final String CASE_ID = "caseID";
        public static final String CULL_DAY = "cullDay";
        public static final String EXAMINATION_DAY = "examinationDay";
        public static final String COORDINATES = "coordinates";
        public static final String INFECTION_TIME_BRANCH_POSITION = "infectionTimeBranchPosition";
        public static final String INFECTIOUS_PERIOD_DISTRIBUTION = "infectiousPeriodDistribution";

        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            final Parameter riemannSampleSize = (Parameter) xo.getElementFirstChild(RIEMANN_SAMPLE_SIZE);
            final boolean hasGeography = xo.hasAttribute(HAS_GEOGRAPHY)
                    ? (Boolean) xo.getAttribute(HAS_GEOGRAPHY) : false;
            final Taxa taxa = (Taxa) xo.getChild(Taxa.class);
            CategoryOutbreak cases = new CategoryOutbreak(null, taxa, hasGeography, riemannSampleSize);
            for(int i=0; i<xo.getChildCount(); i++){
                Object cxo = xo.getChild(i);
                if(cxo instanceof XMLObject && ((XMLObject)cxo).getName().equals(CategoryCase.CATEGORY_CASE)){
                    parseCase((XMLObject)cxo, cases);
                }
            }
            return cases;
        }

        public void parseCase(XMLObject xo, CategoryOutbreak outbreak)
                throws XMLParseException {
            String farmID = (String) xo.getAttribute(CASE_ID);
            final Date cullDate = (Date) xo.getElementFirstChild(CULL_DAY);
            final Date examDate = (Date) xo.getElementFirstChild(EXAMINATION_DAY);
            final ParametricDistributionModel infectiousDist =
                    (ParametricDistributionModel)xo.getElementFirstChild(INFECTIOUS_PERIOD_DISTRIBUTION);
            final Parameter coords = xo.hasChildNamed(COORDINATES) ?
                    (Parameter) xo.getElementFirstChild(COORDINATES) : null;
            Taxa taxa = new Taxa();
            for(int i=0; i<xo.getChildCount(); i++){
                if(xo.getChild(i) instanceof Taxon){
                    taxa.addTaxon((Taxon)xo.getChild(i));
                }
            }
            outbreak.addCase(farmID, examDate, cullDate, infectiousDist, coords, taxa);
        }



        public String getParserDescription(){
            return "Parses a set of 'category' farm cases and the information that they all share";
        }

        public Class getReturnType(){
            return SimpleOutbreak.class;
        }

        public String getParserName(){
            return CATEGORY_OUTBREAK;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] caseRules = {
                new StringAttributeRule(CASE_ID, "The unique identifier for this farm"),
                new ElementRule(CULL_DAY, Date.class, "The date this farm was culled", false),
                new ElementRule(EXAMINATION_DAY, Date.class, "The date this farm was examined", false),
                new ElementRule(Taxon.class, 0, Integer.MAX_VALUE),
                new ElementRule(INFECTIOUS_PERIOD_DISTRIBUTION, ParametricDistributionModel.class, "The probability" +
                        "distribution from which the infectious period of this case is drawn"),
                new ElementRule(INFECTION_TIME_BRANCH_POSITION, Parameter.class, "The exact position on the branch" +
                        " along which the infection of this case occurs that it actually does occur"),
                new ElementRule(COORDINATES, Parameter.class, "The spatial coordinates of this case", true)
        };

        private final XMLSyntaxRule[] rules = {
                new ElementRule(ProductStatistic.class, 0,2),
                new ElementRule(CategoryCase.CATEGORY_CASE, caseRules, 1, Integer.MAX_VALUE),
                new ElementRule(RIEMANN_SAMPLE_SIZE, Parameter.class, "The sample size for the Riemann numerical" +
                        "integration method, used by all child cases.", true),
                new ElementRule(Taxa.class),
                new ElementRule(INFECTIOUS_PERIOD_DISTRIBUTIONS, ParametricDistributionModel.class,
                        "One or more probability distributions for the infectious periods of cases in the oubreak", 1,
                        Integer.MAX_VALUE),
                AttributeRule.newBooleanRule(HAS_GEOGRAPHY, true)
        };
    };


}
