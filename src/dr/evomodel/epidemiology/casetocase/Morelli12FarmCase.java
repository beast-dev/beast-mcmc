package dr.evomodel.epidemiology.casetocase;

import dr.evolution.util.Date;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evolution.util.Units;
import dr.inference.distribution.GammaDistributionModel;
import dr.inference.distribution.ParametricDistributionModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.math.IntegrableUnivariateFunction;
import dr.math.RiemannApproximation;
import dr.math.UnivariateFunction;
import dr.xml.*;

/**
 * Adaptation of the farm incubation and infectious period models from Morelli et al, PLoS Computational Biology, 2012
 * (10.1371/journal.pcbi.1002768.g001)
 *
 * User: Matthew Hall
 * Date: 07/09/2012
 * Time: 16:17
 */

public class Morelli12FarmCase extends AbstractCase {

    public Morelli12FarmCase(String name, String caseID, Date examDate, Date cullDate, Parameter oldestLesionAge,
                             Parameter d, ParametricDistributionModel incubationPeriod, Taxa associatedTaxa){
        super(name);
        this.caseID = caseID;
        //The time value for end of these days is the numerical value of these dates plus 1.
        this.examDate = examDate;
        this.d = d;
        endOfInfectiousDate = cullDate;
        this.associatedTaxa = associatedTaxa;
        this.oldestLesionAge = oldestLesionAge;
        this.incubationPeriod = incubationPeriod;
        infectionDate = new InfectionDatePDF();
        rebuildInfDistribution();
        addModel(infectiousPeriod);
        addModel(this.incubationPeriod);
        addVariable(this.d);

    }


    public Morelli12FarmCase(String caseID, Date examDate, Date cullDate, Parameter oldestLesionAge, Parameter d,
                             ParametricDistributionModel incubationPeriod, Taxa associatedTaxa){
        this(MORELLI_12_FARM_CASE, caseID, examDate, cullDate, oldestLesionAge, d, incubationPeriod,
                associatedTaxa);
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

    public void installNumericalIntegrator(RiemannApproximation integrator) {
        numericalIntegrator = integrator;
    }

    public Taxa getAssociatedTaxa() {
        return associatedTaxa;
    }

    public void setIncubationPeriodDistribution(ParametricDistributionModel distribution){
        incubationPeriod = distribution;
    }

    public ParametricDistributionModel getIncubationPeriodDistribution(){
        return incubationPeriod;
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

    /*Parser.*/

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public static final String CASE_ID = "caseID";
        public static final String CULL_DAY = "cullDay";
        public static final String EXAMINATION_DAY = "examinationDay";
        public static final String INCUBATION_PERIOD_DISTRIBUTION = "incubationPeriodDistribution";
        public static final String OLDEST_LESION_AGE = "oldestLesionAge";
        public static final String SQRT_INFECTIOUS_SCALE = "sqrtInfectiousScale";

        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            String farmID = (String) xo.getAttribute(CASE_ID);
            final Date cullDate = (Date) xo.getElementFirstChild(CULL_DAY);
            final Date examinationDate = (Date) xo.getElementFirstChild(EXAMINATION_DAY);
            final Parameter oldestLesionAgeParameter = (Parameter) xo.getElementFirstChild(OLDEST_LESION_AGE);
            final Parameter d = (Parameter) xo.getElementFirstChild(SQRT_INFECTIOUS_SCALE);
            Taxa tempTaxa = new Taxa();
            for(int i=0; i<xo.getChildCount(); i++){
                if(xo.getChild(i) instanceof Taxon){
                    tempTaxa.addTaxon((Taxon)xo.getChild(i));
                }
            }
            final Taxa associatedTaxa = tempTaxa;
            final ParametricDistributionModel incubationPeriodDistribution =
                    (ParametricDistributionModel) xo.getElementFirstChild(INCUBATION_PERIOD_DISTRIBUTION);
            return new Morelli12FarmCase(farmID, examinationDate, cullDate, oldestLesionAgeParameter, d,
                    incubationPeriodDistribution, associatedTaxa);
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                new StringAttributeRule(CASE_ID, "The unique identifier for this farm"),
                new ElementRule(CULL_DAY, Date.class, "The date this farm was culled", false),
                new ElementRule(EXAMINATION_DAY, Date.class, "The date this farm was examined", false),
                new ElementRule(Taxon.class, 0, Integer.MAX_VALUE),
                new ElementRule(SQRT_INFECTIOUS_SCALE, Parameter.class, "The square root of the scale parameter of " +
                        "all infectiousness periods (variances are proportional to the square of this, see Morelli" +
                        "2012).", false),
                new ElementRule(INCUBATION_PERIOD_DISTRIBUTION, ParametricDistributionModel.class, "The probability " +
                        "distribution of the incubation period of this farm", false),
                new ElementRule(OLDEST_LESION_AGE, Parameter.class, "The estimated oldest lesion date as determined" +
                        "by investigating vets")
        };

        @Override
        public String getParserDescription() {
            return "Parses XML (generated by oldMacDonald.py - contact mdhall272@gmail.com if interested) for " +
                    "descriptions of cases. Incubation periods are gamma distributed, dates of infectiousness have a " +
                    "truncated normal distribution with a maximum value equal to the time of sampling. Infection date" +
                    " is estimated by adding these together.";
        }

        @Override
        public Class getReturnType() {
            return Morelli12FarmCase.class;
        }

        public String getParserName() {
            return MORELLI_12_FARM_CASE;
        }
    };

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

    public class JointDistribution implements IntegrableUnivariateFunction{

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

    public class InfectionDatePDF implements UnivariateFunction{

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
    private Parameter d;
    private ParametricDistributionModel infectiousDate;
    private InfectionDatePDF infectionDate;
    private Taxa associatedTaxa;
    private RiemannApproximation numericalIntegrator;
    private ParametricDistributionModel incubationPeriod;
    private ParametricDistributionModel infectiousPeriod;
    private ParametricDistributionModel storedInfectiousPeriod;
    private double[] fastInfectionDate;

}
