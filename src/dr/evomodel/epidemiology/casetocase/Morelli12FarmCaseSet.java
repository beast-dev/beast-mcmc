package dr.evomodel.epidemiology.casetocase;

import dr.evolution.util.Date;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evolution.util.Units;
import dr.inference.distribution.GammaDistributionModel;
import dr.inference.distribution.ParametricDistributionModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.ProductStatistic;
import dr.inference.model.Variable;
import dr.math.IntegrableUnivariateFunction;
import dr.math.RiemannApproximation;
import dr.math.UnivariateFunction;
import dr.xml.*;

import java.util.ArrayList;

/**
 * Adaptation of the farm incubation and infectious period models from Morelli et al, PLoS Computational Biology, 2012
 * (10.1371/journal.pcbi.1002768.g001)
 *
 * User: Matthew Hall
 * Date: 07/09/2012
 * Time: 16:17
 */
public class Morelli12FarmCaseSet extends AbstractCaseSet{

    public Morelli12FarmCaseSet(String name, ParametricDistributionModel incubationPeriodDistribution, Parameter d,
                                ArrayList<AbstractCase> farms, Parameter riemannSampleSize){
        this(name,incubationPeriodDistribution,d,riemannSampleSize);
        this.cases = farms;
        for(AbstractCase farm : farms){
            addModel(farm);
        }
    }

    public Morelli12FarmCaseSet(ParametricDistributionModel incubationPeriodDistribution, Parameter d,
                                ArrayList<AbstractCase> farms, Parameter riemannSampleSize){
        this(MORELLI_12_FARM_CASE_SET, incubationPeriodDistribution, d, farms, riemannSampleSize);
    }

    // with the inner class, initialisation has to take places without cases - add them later

    public Morelli12FarmCaseSet(String name, ParametricDistributionModel incubationPeriodDistribution, Parameter d,
                                Parameter riemannSampleSize){
        super(name);
        this.incubationPeriod = incubationPeriodDistribution;
        addModel(this.incubationPeriod);
        this.d = d;
        numericalIntegrator = new RiemannApproximation((int)riemannSampleSize.getParameterValue(0));
        cases = new ArrayList<AbstractCase>();
    }

    public Morelli12FarmCaseSet(ParametricDistributionModel incubationPeriodDistribution, Parameter d,
                                Parameter riemannSampleSize){
        this(MORELLI_12_FARM_CASE_SET, incubationPeriodDistribution, d, riemannSampleSize);
    }

    private void addCase(String caseID, Date examDate, Date cullDate, Parameter oldestLesionAge, Taxa associatedTaxa){
        Morelli12FarmCase thisCase = new Morelli12FarmCase(caseID, examDate, cullDate, oldestLesionAge, associatedTaxa);
        cases.add(thisCase);
    }

    /* Likelihood of the root branch (the farm is infectious by the root node time).*/

    public double noTransmissionBranchLikelihood(AbstractCase farm, Integer farmInfectiousBy) {
        return Math.exp(noTransmissionBranchLogLikelihood(farm, farmInfectiousBy));
    }

    public double noTransmissionBranchLogLikelihood(AbstractCase farm, Integer farmInfectiousBy) {
        if(farm.culledYet(farmInfectiousBy)){
            return Double.NEGATIVE_INFINITY;
        } else {
            return Math.log(((Morelli12FarmCase) farm).infectiousCDF(farmInfectiousBy));
        }
    }

    /* Likelihood of a non-root branch (the farm is infected at the parent node time and infectious by the child node
    time). */

    public double transmissionBranchLikelihood(AbstractCase parent, AbstractCase child, Integer childInfected, Integer
            childInfectiousBy) {
        return Math.exp(transmissionBranchLogLikelihood(parent, child, childInfected, childInfectiousBy));
    }

    public double transmissionBranchLogLikelihood(AbstractCase parent, AbstractCase child, Integer childInfected, Integer
            childInfectiousBy) {
        if(child.culledYet(childInfectiousBy)){
            return Double.NEGATIVE_INFINITY;
        } else if(parent==child){
            return 0;
        } else {
            return Math.log(((Morelli12FarmCase)child).periodInfectionDistribution(childInfected - 1, childInfected,
                    childInfectiousBy));
        }
    }

    public ArrayList<AbstractCase> getCases() {
        return new ArrayList<AbstractCase>(cases);
    }

    /* Parser. */

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {
        //for the case set

        public static final String INCUBATION_PERIOD_DISTRIBUTION = "incubationPeriod";
        public static final String RIEMANN_SAMPLE_SIZE = "riemannSampleSize";
        public static final String SQRT_INFECTIOUS_SCALE = "sqrtInfectiousScale";
        public static final String CASE_NAME = "morelli12FarmCase";

        //for the cases

        public static final String CASE_ID = "caseID";
        public static final String CULL_DAY = "cullDay";
        public static final String EXAMINATION_DAY = "examinationDay";
        public static final String OLDEST_LESION_AGE = "oldestLesionAge";

        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            final ParametricDistributionModel incubationPeriodDistribution =
                    (ParametricDistributionModel) xo.getElementFirstChild(INCUBATION_PERIOD_DISTRIBUTION);
            final Parameter d = (Parameter) xo.getElementFirstChild(SQRT_INFECTIOUS_SCALE);
            final Parameter riemannSampleSize = (Parameter) xo.getElementFirstChild(RIEMANN_SAMPLE_SIZE);
            Morelli12FarmCaseSet cases = new Morelli12FarmCaseSet(incubationPeriodDistribution, d, riemannSampleSize);
            for(int i=0; i<xo.getChildCount(); i++){
                if(xo.getName().equals(CASE_NAME)){
                    parseCase((XMLObject)xo.getChild(i),cases);
                }
            }
            return cases;
        }

        public void parseCase(XMLObject xo, Morelli12FarmCaseSet caseSet) throws XMLParseException {
            String farmID = (String) xo.getAttribute(CASE_ID);
            final Date cullDate = (Date) xo.getElementFirstChild(CULL_DAY);
            final Date examDate = (Date) xo.getElementFirstChild(EXAMINATION_DAY);
            final Parameter oldestLesionAge = (Parameter) xo.getElementFirstChild(OLDEST_LESION_AGE);
            Taxa taxa = new Taxa();
            for(int i=0; i<xo.getChildCount(); i++){
                if(xo.getChild(i) instanceof Taxon){
                    taxa.addTaxon((Taxon)xo.getChild(i));
                }
            }
            caseSet.addCase(farmID, examDate, cullDate, oldestLesionAge, taxa);
        }

        @Override
        public String getParserDescription(){
            return "Parses a set of Morelli 2012 farm cases and the information that they all share";
        }

        @Override
        public Class getReturnType(){
            return Morelli12FarmCaseSet.class;
        }

        public String getParserName(){
            return MORELLI_12_FARM_CASE_SET;
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] caseRules = {
                new StringAttributeRule(CASE_ID, "The unique identifier for this farm"),
                new ElementRule(CULL_DAY, Date.class, "The date this farm was culled", false),
                new ElementRule(EXAMINATION_DAY, Date.class, "The date this farm was examined", false),
                new ElementRule(Taxon.class, 0, Integer.MAX_VALUE),
                new ElementRule(OLDEST_LESION_AGE, Parameter.class, "The estimated oldest lesion date as determined" +
                        "by investigating vets")
        };

        private final XMLSyntaxRule[] rules = {
                new ElementRule(ProductStatistic.class, 0,2),
                new ElementRule(INCUBATION_PERIOD_DISTRIBUTION, ParametricDistributionModel.class, "The probability " +
                        "distribution of incubation periods (constructed in the XML so farm elements can inherit" +
                        "it).", false),
                new ElementRule(CASE_NAME, caseRules, 1, Integer.MAX_VALUE),
                new ElementRule(SQRT_INFECTIOUS_SCALE, Parameter.class, "The square root of the scale parameter of " +
                        "all infectiousness periods (variances are proportional to the square of this, see Morelli" +
                        "2012).", false),
                new ElementRule(RIEMANN_SAMPLE_SIZE, Parameter.class, "The sample size for the Riemann numerical" +
                        "integration method, used by all child cases.", true)
        };
    };

    public static final String MORELLI_12_FARM_CASE_SET = "morelli12FarmCaseSet";
    private ParametricDistributionModel incubationPeriod;
    private RiemannApproximation numericalIntegrator;
    private Parameter d;

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        fireModelChanged();
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
        fireModelChanged();
    }

    @Override
    protected void storeState() {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    protected void restoreState() {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    protected void acceptState() {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    //Case class.

    private class Morelli12FarmCase extends AbstractCase {
        public Morelli12FarmCase(String name, String caseID, Date examDate, Date cullDate, Parameter oldestLesionAge,
                                           Taxa associatedTaxa){
            super(name);
            this.caseID = caseID;
            //The time value for end of these days is the numerical value of these dates plus 1.
            this.examDate = examDate;
            endOfInfectiousDate = cullDate;
            this.associatedTaxa = associatedTaxa;
            this.oldestLesionAge = oldestLesionAge;
            infectionDate = new InfectionDatePDF();
            rebuildInfDistribution();
            addModel(infectiousPeriod);
            addModel(incubationPeriod);
            addVariable(d);

        }


        public Morelli12FarmCase(String caseID, Date examDate, Date cullDate, Parameter oldestLesionAge,
                                 Taxa associatedTaxa){
            this(MORELLI_12_FARM_CASE, caseID, examDate, cullDate, oldestLesionAge, associatedTaxa);
        }

        private void rebuildInfDistribution(){
            Parameter infectious_shape = new Parameter.Default
                    (oldestLesionAge.getParameterValue(0)/Math.pow(d.getParameterValue(0),2));
            Parameter infectious_scale = new Parameter.Default(Math.pow(d.getParameterValue(0),2));
            infectiousPeriod = new GammaDistributionModel(infectious_shape, infectious_scale);
        }

        public Date getLatestPossibleInfectionDate() {
            Double doubleDate = examDate.getTimeValue();
            return Date.createTimeSinceOrigin(doubleDate, Units.Type.DAYS, examDate.getOrigin());
        }

        public Taxa getAssociatedTaxa() {
            return associatedTaxa;
        }

        public void setInfectiousPeriodDistribution(ParametricDistributionModel distribution){
            infectiousPeriod = distribution;
        }

        public ParametricDistributionModel getInfectiousPeriodDistribution(){
            return infectiousPeriod;
        }

        public boolean culledYet(int day) {
            return day>endOfInfectiousDate.getTimeValue()+1;
        }

        public Object getInfectionDate() {
            return infectionDate;
        }

        public Object getInfectiousDate() {
            return null;
        }

        public Object getEndOfInfectiousDate() {
            return endOfInfectiousDate;
        }

        public Double getEndOfInfectiousDateModeHeight(Date latestTaxonDate) {
            return latestTaxonDate.getTimeValue()-endOfInfectiousDate.getTimeValue();
        }

        public Double getInfectiousDateModeHeight(Date latestTaxonDate) {
            return null;
        }

        public double infectiousCDF(double time){
            return 1-infectiousPeriod.cdf(oldestLesionAge.getParameterValue(0)-time);
        }

    /* Probability that infection occurred between 'earliestInfection' and 'latestInfection' given that the case is
    infectious by 'infectiousBy' (NOT that it actually becomes infectious at that time). */

        public double periodInfectionDistribution(double earliestInfection, double latestInfection, double infectiousBy){
            JointDistribution tempDist = new JointDistribution(infectiousBy);
            return tempDist.evaluateIntegral(earliestInfection, latestInfection);
        }

    /* Probability distribution of infection date being 'argument' given that the case is infectious by
    'infectiousBy'. (NOT that it actually becomes infectious at that time). */

        public double infectionDistribution(double argument, double infectiousBy){
            JointDistribution tempDist = new JointDistribution(infectiousBy);
            return tempDist.evaluate(argument);
        }

        @Override
        protected void handleModelChangedEvent(Model model, Object object, int index) {
            rebuildInfDistribution();
            fireModelChanged();
        }

        @Override
        protected void handleVariableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
            rebuildInfDistribution();
            fireModelChanged();
        }

        @Override
        protected void storeState() {
            storedInfectiousPeriod = infectiousPeriod;
        }

        @Override
        protected void restoreState() {
            infectiousPeriod = storedInfectiousPeriod;
        }

        @Override
        protected void acceptState() {
            //To change body of implemented methods use File | Settings | File Templates.
        }

    /* Probability for the infection date taking the value 'argument' and the incubation period being at most
     * 't2' minus 'argument'; i.e. infected at 'argument' and infectious before 't2'. */

        public class JointDistribution implements IntegrableUnivariateFunction {

            public JointDistribution(double t2){
                this.t2=t2;
            }

            public double evaluateIntegral(double a, double b) {
                return numericalIntegrator.integrate(this, a, b);
            }

            public double evaluate(double argument) {
                return infectionDate.evaluate(argument)*incubationPeriod.cdf(t2-argument);
            }

            public double getLowerBound() {
                return Double.NEGATIVE_INFINITY;
            }

            public double getUpperBound() {
                return t2;
            }

            private double t2;
        }


    /* Inner class for the probability distribution of the sum of the infectious date and minus the incubation period,
    i.e. the infection date. Not implementing Distribution since calculation of the mean, etc. is unnecessary.
     */

        public class InfectionDatePDF implements UnivariateFunction {

            public InfectionDatePDF(){
            }

            public double evaluate(double argument) {
                if(numericalIntegrator==null){
                    throw new RuntimeException("Numerical integrator not specified.");
                } else {
                    L_x_fi underlyingFunction = new L_x_fi(argument);
                    return numericalIntegrator.integrate(underlyingFunction,0,examDate.getTimeValue()-argument);
                }
            }

            public double getLowerBound() {
                return Double.NEGATIVE_INFINITY;
            }

            public double getUpperBound() {
                return Double.POSITIVE_INFINITY;
            }
        }

    /* Calculation of InfectionDatePDF requires another integral */

        public class L_x_fi implements UnivariateFunction{

            public L_x_fi(double currentT){
                this.currentT = currentT;
            }

            public double evaluate(double argument) {
                if(argument<0 || argument>examDate.getTimeValue() - currentT){
                    return 0;
                } else {
                    return incubationPeriod.pdf(examDate.getTimeValue() - currentT - argument)
                            *infectiousPeriod.pdf(argument);
                }
            }

            public double getLowerBound() {
                return Double.NEGATIVE_INFINITY;
            }

            public double getUpperBound() {
                return Double.POSITIVE_INFINITY;
            }
            private double currentT;
        }

        public static final String MORELLI_12_FARM_CASE = "morelli12FarmCase";
        private Date examDate;
        private Date endOfInfectiousDate;
        private Parameter oldestLesionAge;
        private ParametricDistributionModel infectiousDate;
        private InfectionDatePDF infectionDate;
        private Taxa associatedTaxa;
        private ParametricDistributionModel infectiousPeriod;
        private ParametricDistributionModel storedInfectiousPeriod;
    }

}
