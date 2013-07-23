package dr.evomodel.epidemiology.casetocase;

import dr.evolution.tree.Tree;
import dr.inference.model.*;
import dr.math.IntegrableUnivariateFunction;
import dr.math.Integral;
import dr.math.RiemannApproximation;
import dr.math.UnivariateFunction;
import dr.xml.*;

import java.util.*;

/**
 * A likelihood function for transmission between identified epidemiological cases
 *
 * Timescale must be in days. Python scripts to write XML for it and analyse the posterior set of networks exist;
 * contact MH. @todo make timescale not just in days
 *
 * Latent periods are not implemented currently
 *
 * @author Matthew Hall
 * @version $Id: $
 */

public class CaseToCaseTransmissionLikelihood extends AbstractModelLikelihood {

    private AbstractOutbreak outbreak;
    private CaseToCaseTreeLikelihood treeLikelihood;
    private SpatialKernel spatialKernel;
    private Parameter kernelAlpha;
    private Parameter transmissionRate;
    private boolean hasLatentPeriods;
    private Integer[] parents;
    private Integer[] storedParents;
    private HashMap<Integer, Double> infectionTimes;
    private HashMap<Integer, Double> storedInfectionTimes;
    private InnerIntegral probFunct;
    private boolean likelihoodKnown;
    private double logLikelihood;
    public static final String CASE_TO_CASE_TRANSMISSION_LIKELIHOOD = "caseToCaseTransmissionLikelihood";

    public CaseToCaseTransmissionLikelihood(String name, AbstractOutbreak outbreak,
                                            CaseToCaseTreeLikelihood treeLikelihood, SpatialKernel spatialKernal,
                                            Parameter transmissionRate){
        super(name);
        this.outbreak = outbreak;
        this.treeLikelihood = treeLikelihood;
        this.spatialKernel = spatialKernal;
        kernelAlpha = spatialKernal.geta();
        this.transmissionRate = transmissionRate;
        hasLatentPeriods = outbreak.hasLatentPeriods();
        this.addModel(treeLikelihood);
        this.addVariable(transmissionRate);
        parents = treeLikelihood.getParentsArray();
        infectionTimes = treeLikelihood.getInfTimesMap();
        probFunct = new InnerIntegral();
        likelihoodKnown = false;
    }


    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if(model instanceof CaseToCaseTreeLikelihood){
            parents = treeLikelihood.getParentsArray();
            infectionTimes = treeLikelihood.getInfTimesMap();
            Arrays.sort(parents, new CaseNoInfectionComparator());
        }
        likelihoodKnown = false;
        fireModelChanged();
    }

    @Override
    protected void storeState() {
        storedParents = Arrays.copyOf(parents, parents.length);
        storedInfectionTimes = new HashMap<Integer, Double>(infectionTimes);
    }

    @Override
    protected void restoreState() {
        parents = storedParents;
        infectionTimes = storedInfectionTimes;
    }

    @Override
    protected void acceptState() {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        likelihoodKnown = false;
        fireModelChanged();
    }

    public Model getModel() {
        return this;
    }

    public double getLogLikelihood() {
        if(!likelihoodKnown){
            int N = outbreak.size();
            double lambda = transmissionRate.getParameterValue(0);
            double unnormalisedValue = Math.pow(lambda, N-1)*aAlpha()*Math.exp(-lambda*bAlpha())/N;
            double normalisationValue = probFunct.evaluateIntegral(0, kernelAlpha.getBounds().getUpperLimit(0));
            double treeLogL = treeLikelihood.getLogLikelihood();
            likelihoodKnown = true;
            logLikelihood =  Math.log(treeLogL) + Math.log(unnormalisedValue) - Math.log(normalisationValue);
        }
        return logLikelihood;
    }


    public void makeDirty() {
        likelihoodKnown = false;
    }

    private class CaseInfectionComparator implements Comparator<AbstractCase> {
        public int compare(AbstractCase abstractCase, AbstractCase abstractCase2) {
            return Double.compare(infectionTimes.get(outbreak.getCaseIndex(abstractCase)),
                    infectionTimes.get(outbreak.getCaseIndex(abstractCase2)));
        }
    }

    private class CaseNoInfectionComparator implements Comparator<Integer> {
        public int compare(Integer integer, Integer integer2) {
            return Double.compare(infectionTimes.get(integer), infectionTimes.get(integer2));
        }
    }

    private double bAlpha(){
        double total = 0;
        ArrayList<AbstractCase> copyOfCases = new ArrayList<AbstractCase>(outbreak.getCases());
        Collections.sort(copyOfCases, new CaseInfectionComparator());
        for(int i=1; i<outbreak.size(); i++){
            for(int j=0; j<i; j++){
                total += (infectionTimes.get(outbreak.getCaseIndex(copyOfCases.get(i))) -
                        infectionTimes.get(outbreak.getCaseIndex(copyOfCases.get(j))))
                        * outbreak.getKernalValue(copyOfCases.get(i), copyOfCases.get(j), spatialKernel);
            }
        }
        return total;
    }

    private double bAlpha(double alpha){
        double total = 0;
        ArrayList<AbstractCase> copyOfCases = new ArrayList<AbstractCase>(outbreak.getCases());
        Collections.sort(copyOfCases, new CaseInfectionComparator());
        for(int i=1; i<outbreak.size(); i++){
            for(int j=0; j<i; j++){
                total += (infectionTimes.get(outbreak.getCaseIndex(copyOfCases.get(i))) -
                        infectionTimes.get(outbreak.getCaseIndex(copyOfCases.get(j))))
                        * outbreak.getKernalValue(copyOfCases.get(i), copyOfCases.get(j), spatialKernel, alpha);
            }
        }
        return total;
    }

    private double aAlpha(){
        double product = 1;
        ArrayList<AbstractCase> copyOfCases = new ArrayList<AbstractCase>(outbreak.getCases());
        Collections.sort(copyOfCases, new CaseInfectionComparator());
        for(int i=1; i<outbreak.size(); i++){
            product *= outbreak.getKernalValue(copyOfCases.get(i), copyOfCases.get(i=1), spatialKernel);
        }
        return product;
    }

    private double aAlpha(double alpha){
        double product = 1;
        ArrayList<AbstractCase> copyOfCases = new ArrayList<AbstractCase>(outbreak.getCases());
        Collections.sort(copyOfCases, new CaseInfectionComparator());
        for(int i=1; i<outbreak.size(); i++){
            product *= outbreak.getKernalValue(copyOfCases.get(i), copyOfCases.get(i=1), spatialKernel, alpha);
        }
        return product;
    }

    // the integral in terms of only alpha

    private class InnerIntegral implements IntegrableUnivariateFunction {

        RiemannApproximation finiteIntegrator;
        InfiniteUpperLimitRiemannApproximation infiniteIntegrator;

        private InnerIntegral(){
            finiteIntegrator = new RiemannApproximation(50);
            infiniteIntegrator = new InfiniteUpperLimitRiemannApproximation(1, 0.001);
        }


        @Override
        public double evaluateIntegral(double a, double b) {
            if(b!=Double.POSITIVE_INFINITY){
                return finiteIntegrator.integrate(this, a, b);
            } else {
                return infiniteIntegrator.integrate(this, a, b);
            }
        }

        @Override
        public double evaluate(double argument) {
            return aAlpha(argument)/Math.pow(bAlpha(argument), outbreak.size());
        }

        @Override
        public double getLowerBound() {
            return 0;
        }

        @Override
        public double getUpperBound() {
            return Double.POSITIVE_INFINITY;
        }
    }

    //rudimentary numerical integrator for infinite limit. Presumes the limit of the function is 0.

    private class InfiniteUpperLimitRiemannApproximation implements Integral {

        final double step;
        final double tolerance;

        private InfiniteUpperLimitRiemannApproximation(double step, double tolerance){
            this.step = step;
            this.tolerance = tolerance;
        }

        public double integrate(UnivariateFunction f, double min, double max) {
            if(max!=Double.POSITIVE_INFINITY){
                throw new RuntimeException("Shouldn't be using this integrator for this integral");
            }
            double integral = 0;
            double gridpoint = min;
            boolean farEnough = false;
            while(!farEnough){
                double increase = f.evaluate(gridpoint);
                integral += increase*step;
                gridpoint += step;
                if(Math.abs(increase)<tolerance){
                    farEnough = true;
                }
            }
            return integral;
        }
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public static final String TRANSMISSION_RATE = "transmissionRate";

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            CaseToCaseTreeLikelihood c2cTL = (CaseToCaseTreeLikelihood)
                    xo.getChild(CaseToCaseTreeLikelihood.class);
            SpatialKernel kernel = (SpatialKernel) xo.getChild(SpatialKernel.class);
            Parameter transmissionRate = (Parameter) xo.getElementFirstChild(TRANSMISSION_RATE);
            return new CaseToCaseTransmissionLikelihood(CASE_TO_CASE_TRANSMISSION_LIKELIHOOD, c2cTL.getOutbreak(),
                    c2cTL, kernel, transmissionRate);
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        public String getParserDescription() {
            return "This element represents a probability distribution for epidemiological parameters of an outbreak" +
                    "given a phylogenetic tree";
        }

        public Class getReturnType() {
            return CaseToCaseTransmissionLikelihood.class;
        }

        public String getParserName() {
            return CASE_TO_CASE_TRANSMISSION_LIKELIHOOD;
        }

        private final XMLSyntaxRule[] rules = {
                new ElementRule(CaseToCaseTreeLikelihood.class, "The tree likelihood"),
                new ElementRule(SpatialKernel.class, "The spatial kernel"),
                new ElementRule(Parameter.class, "The transmission rate")
        };

    };



}


