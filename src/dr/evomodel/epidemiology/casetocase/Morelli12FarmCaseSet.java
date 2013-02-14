package dr.evomodel.epidemiology.casetocase;

import dr.inference.distribution.ParametricDistributionModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.ProductStatistic;
import dr.inference.model.Variable;
import dr.math.RiemannApproximation;
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
        super(name);

        this.incubationPeriodDistribution = incubationPeriodDistribution;
        addModel(this.incubationPeriodDistribution);
        this.cases = farms;
        this.d = d;
        numericalIntegrator = new RiemannApproximation((int)riemannSampleSize.getParameterValue(0));
        for(AbstractCase farm : farms){
            ((Morelli12FarmCase)farm).installNumericalIntegrator(numericalIntegrator);
            addModel(farm);
        }
    }

    public Morelli12FarmCaseSet(ParametricDistributionModel incubationPeriodDistribution, Parameter d,
                                ArrayList<AbstractCase> farms, Parameter riemannSampleSize){
        this(MORELLI_12_FARM_CASE_SET, incubationPeriodDistribution, d, farms, riemannSampleSize);
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
        //The duplication, believe it or not, makes sense since the distribution needs to exist in the XML to
        // let the farm elements inherit it, but the operators should apply to the mean and stdev.
        public static final String INCUBATION_PERIOD_MEAN = "incubationPeriodMean";
        public static final String INCUBATION_PERIOD_STDEV = "incubationPeriodStdev";
        public static final String INCUBATION_PERIOD_DISTRIBUTION = "incubationPeriodDistribution";
        public static final String RIEMANN_SAMPLE_SIZE = "riemannSampleSize";
        public static final String SQRT_INFECTIOUS_SCALE = "sqrtInfectiousScale";

        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            final ParametricDistributionModel incubationPeriodDistribution =
                    (ParametricDistributionModel) xo.getElementFirstChild(INCUBATION_PERIOD_DISTRIBUTION);
            ArrayList<AbstractCase> tempFarms = new ArrayList<AbstractCase>();
            for(int i=0; i<xo.getChildCount(); i++){
                if(xo.getChild(i) instanceof Morelli12FarmCase){
                    tempFarms.add((Morelli12FarmCase)xo.getChild(i));
                }
            }
            final Parameter d = (Parameter) xo.getElementFirstChild(SQRT_INFECTIOUS_SCALE);
            final Parameter riemannSampleSize = (Parameter) xo.getElementFirstChild(RIEMANN_SAMPLE_SIZE);
            final ArrayList<AbstractCase> farms = tempFarms;
            return new Morelli12FarmCaseSet(incubationPeriodDistribution, d, farms, riemannSampleSize);
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

        private final XMLSyntaxRule[] rules = {
                new ElementRule(ProductStatistic.class, 0,2),
                new ElementRule(INCUBATION_PERIOD_STDEV, ProductStatistic.class, "The standard deviation of the distribution" +
                        "of incubation periods", true),
                new ElementRule(INCUBATION_PERIOD_DISTRIBUTION, ParametricDistributionModel.class, "The probability " +
                        "distribution of incubation periods (constructed in the XML so farm elements can inherit" +
                        "it).", false),
                new ElementRule(Morelli12FarmCase.class, 1, Integer.MAX_VALUE),
                new ElementRule(SQRT_INFECTIOUS_SCALE, Parameter.class, "The square root of the scale parameter of " +
                        "all infectiousness periods (variances are proportional to the square of this, see Morelli" +
                        "2012).", false),
                new ElementRule(RIEMANN_SAMPLE_SIZE, Parameter.class, "The sample size for the Riemann numerical" +
                        "integration method, used by all child cases.", true)
        };
    };

    public static final String MORELLI_12_FARM_CASE_SET = "morelli12FarmCaseSet";
    private ParametricDistributionModel incubationPeriodDistribution;
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
}
