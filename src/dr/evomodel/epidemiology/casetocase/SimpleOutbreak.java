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
import dr.xml.*;

import java.util.ArrayList;

/**
 * Simple version with no latent period; built for speed and extended paintings
 *
 * User: Matthew Hall
 * Date: 07/09/2012
 */

public class SimpleOutbreak extends AbstractOutbreak {

    public SimpleOutbreak(String name, Parameter d, Parameter riemannSampleSize, ArrayList<AbstractCase> farms){
        this(name, d, riemannSampleSize);
        cases = farms;
        for(AbstractCase farm : farms){
            addModel(farm);
        }
    }

    public SimpleOutbreak(Parameter d, ArrayList<AbstractCase> farms, Parameter riemannSampleSize){
        this(SIMPLE_OUTBREAK, d, riemannSampleSize, farms);
    }

    // with the inner class, initialisation has to take places without cases - add them later

    public SimpleOutbreak(String name, Parameter d,
                          Parameter riemannSampleSize){
        super(name);
        this.d = d;
        numericalIntegrator = new RiemannApproximation((int)riemannSampleSize.getParameterValue(0));
        cases = new ArrayList<AbstractCase>();
    }

    public SimpleOutbreak(Parameter d,
                          Parameter riemannSampleSize){
        this(SIMPLE_OUTBREAK, d, riemannSampleSize);
    }

    private void addCase(String caseID, Date examDate, Date cullDate, Parameter oldestLesionAge, Taxa associatedTaxa){
        SimpleCase thisCase = new SimpleCase(caseID, examDate, cullDate, oldestLesionAge, associatedTaxa);
        cases.add(thisCase);
        addModel(thisCase);
    }


    //indices of paintings:
    //0 - node
    //1 - tree.getChild(node,0)
    //2 - tree.getChild(node,1)
    //
    // Similar for times
    // this ignores whether the tree will actually accept the paintings at this node (although it must be a
    // non-extended painting) and just calculates the probability of these timings

    private double localLogP(SimpleCase parentPainting, SimpleCase childPainting, double time1,
                            double time2, boolean extended){
        if(!extended){
            if(parentPainting==childPainting){
                return 0;
            } else {
                //the event in question is that the case corresponding to the child whose painting is different infected
                // at the time of node and infectious by the time of that child.
                return Math.log(childPainting.infectedAt(time1));
            }
        } else {
            if(parentPainting==childPainting){
                return 0;
            } else {
                return Math.log(childPainting.infectedBetween(time1, time2));
            }
        }
    }

    private double localP(SimpleCase parentPainting, SimpleCase childPainting, double time1,
                         double time2, boolean extended){
        return Math.exp(logP(parentPainting, childPainting, time1, time2, extended));
    }

    public double P(AbstractCase parentPainting, AbstractCase childPainting, double time1,
                    double time2, boolean extended){
        return localP((SimpleCase)parentPainting, (SimpleCase)childPainting, time1, time2, extended);
    }

    public double logP(AbstractCase parentPainting, AbstractCase childPainting, double time1,
                       double time2, boolean extended){
        return localLogP((SimpleCase)parentPainting, (SimpleCase)childPainting, time1, time2, extended);
    }

    private double localProbInfectiousBy(SimpleCase painting, double time, boolean extended){
        return Math.exp(logProbInfectiousBy(painting, time, extended));
    }

    private double localLogProbInfectiousBy(SimpleCase painting, double time, boolean extended){
        return Math.log(1-painting.getInfectiousPeriodDistribution().cdf(painting.getEndOfInfectiousPeriod() - time));
    }

    public double probInfectiousBy(AbstractCase painting, double time, boolean extended){
        return localProbInfectiousBy((SimpleCase)painting, time, extended);
    }

    public double logProbInfectiousBy(AbstractCase painting, double time, boolean extended){
        return localLogProbInfectiousBy((SimpleCase)painting, time, extended);
    }

    /* Parser. */

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        //for the outbreak

        public static final String RIEMANN_SAMPLE_SIZE = "riemannSampleSize";
        public static final String SQRT_INFECTIOUS_SCALE = "sqrtInfectiousScale";

        //for the cases

        public static final String CASE_ID = "caseID";
        public static final String CULL_DAY = "cullDay";
        public static final String EXAMINATION_DAY = "examinationDay";
        public static final String OLDEST_LESION_AGE = "oldestLesionAge";

        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            final Parameter d = (Parameter) xo.getElementFirstChild(SQRT_INFECTIOUS_SCALE);
            final Parameter riemannSampleSize = (Parameter) xo.getElementFirstChild(RIEMANN_SAMPLE_SIZE);
            SimpleOutbreak cases = new SimpleOutbreak(d, riemannSampleSize);
            for(int i=0; i<xo.getChildCount(); i++){
                Object cxo = xo.getChild(i);
                if(cxo instanceof XMLObject && ((XMLObject)cxo).getName().equals(SimpleCase.SIMPLE_CASE)){
                    parseCase((XMLObject)cxo, cases);
                }
            }
            return cases;
        }

        public void parseCase(XMLObject xo, SimpleOutbreak outbreak)
                throws XMLParseException {
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
            outbreak.addCase(farmID, examDate, cullDate, oldestLesionAge, taxa);
        }

        @Override
        public String getParserDescription(){
            return "Parses a set of 'simple' farm cases and the information that they all share";
        }

        @Override
        public Class getReturnType(){
            return SimpleOutbreak.class;
        }

        public String getParserName(){
            return SIMPLE_OUTBREAK;
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
                        "by investigating vets"),
        };

        private final XMLSyntaxRule[] rules = {
                new ElementRule(ProductStatistic.class, 0,2),
                new ElementRule(SimpleCase.SIMPLE_CASE, caseRules, 1, Integer.MAX_VALUE),
                new ElementRule(SQRT_INFECTIOUS_SCALE, Parameter.class, "The square root of the scale parameter of " +
                        "all infectiousness periods (variances are proportional to the square of this, see Morelli" +
                        "2012).", false),
                new ElementRule(RIEMANN_SAMPLE_SIZE, Parameter.class, "The sample size for the Riemann numerical" +
                        "integration method, used by all child cases.", true),
        };
    };

    public static final String SIMPLE_OUTBREAK = "simpleOutbreak";
    private final RiemannApproximation numericalIntegrator;
    private final Parameter d;

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

    private class SimpleCase extends AbstractCase {

        // geography version

        public SimpleCase(String name, String caseID, Date examDate, Date cullDate, Parameter oldestLesionAge,
                          Taxa associatedTaxa){
            super(name);
            this.caseID = caseID;
            //The time value for end of these days is the numerical value of these dates plus 1.
            this.examDate = examDate;
            endOfInfectiousDate = cullDate;
            this.associatedTaxa = associatedTaxa;
            this.oldestLesionAge = oldestLesionAge;
            rebuildInfDistribution();
            this.addModel(infectiousPeriodDistribution);
            this.addVariable(d);
        }

        // non-geography version


        public SimpleCase(String caseID, Date examDate, Date cullDate, Parameter oldestLesionAge,
                          Taxa associatedTaxa){
            this(SIMPLE_CASE, caseID, examDate, cullDate, oldestLesionAge, associatedTaxa);
        }


        private void rebuildInfDistribution(){
            Parameter infectious_shape = new Parameter.Default
                    (oldestLesionAge.getParameterValue(0)/Math.pow(d.getParameterValue(0),2));
            Parameter infectious_scale = new Parameter.Default(Math.pow(d.getParameterValue(0),2));
            infectiousPeriodDistribution = new GammaDistributionModel(infectious_shape, infectious_scale);
        }

        public Date getLatestPossibleInfectionDate() {
            Double doubleDate = examDate.getTimeValue();
            return Date.createTimeSinceOrigin(doubleDate, Units.Type.DAYS, examDate.getOrigin());
        }

        public Taxa getAssociatedTaxa() {
            return associatedTaxa;
        }

        public double getEndOfInfectiousPeriod(){
            return endOfInfectiousDate.getTimeValue();
        }

        public void setInfectiousPeriodDistribution(ParametricDistributionModel distribution){
            infectiousPeriodDistribution = distribution;
        }

        public ParametricDistributionModel getInfectiousPeriodDistribution(){
            return infectiousPeriodDistribution;
        }

        public boolean culledYet(int day) {
            return day>endOfInfectiousDate.getTimeValue()+1;
        }

        public boolean culledYet(double time) {
            return time>endOfInfectiousDate.getTimeValue();
        }

        public Object getEndOfInfectiousDate() {
            return endOfInfectiousDate;
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
                double tempEndCDF = infectiousPeriodDistribution.cdf(endOfInfectiousDate.getTimeValue()-endPoint);
                double tempStartCDF = infectiousPeriodDistribution.cdf(endOfInfectiousDate.getTimeValue()-start);

                return infectiousPeriodDistribution.cdf(endOfInfectiousDate.getTimeValue()-start)
                        - infectiousPeriodDistribution.cdf(endOfInfectiousDate.getTimeValue()-endPoint);
            }

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
            storedInfectiousPeriodDistribution = infectiousPeriodDistribution;
        }

        @Override
        protected void restoreState() {
            infectiousPeriodDistribution = storedInfectiousPeriodDistribution;
        }

        @Override
        protected void acceptState() {
        }

        public static final String SIMPLE_CASE = "simpleCase";
        private final Date examDate;
        private final Date endOfInfectiousDate;
        private final Parameter oldestLesionAge;
        private ParametricDistributionModel infectiousPeriodDistribution;
        private ParametricDistributionModel storedInfectiousPeriodDistribution;
    }

}
