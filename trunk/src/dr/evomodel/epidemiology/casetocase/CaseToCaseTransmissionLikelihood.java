package dr.evomodel.epidemiology.casetocase;

import dr.inference.model.*;
import dr.math.IntegrableUnivariateFunction;
import dr.math.Integral;
import dr.math.RiemannApproximation;
import dr.math.UnivariateFunction;
import dr.xml.*;
import org.apache.commons.math.util.MathUtils;

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
    private Integer[] storedParents;
    private HashMap<Integer, Double> infectionTimes;
    private HashMap<Integer, Double> storedInfectionTimes;
    private InnerIntegral probFunct;
    private boolean likelihoodKnown;
    private boolean storedLikelihoodKnown;
    private double logLikelihood;
    private double storedLogLikelihood;
    private boolean hasGeography;
    public static final String CASE_TO_CASE_TRANSMISSION_LIKELIHOOD = "caseToCaseTransmissionLikelihood";

    public CaseToCaseTransmissionLikelihood(String name, AbstractOutbreak outbreak,
                                            CaseToCaseTreeLikelihood treeLikelihood, SpatialKernel spatialKernal,
                                            Parameter transmissionRate){
        super(name);
        this.outbreak = outbreak;
        this.treeLikelihood = treeLikelihood;
        this.spatialKernel = spatialKernal;
        if(spatialKernal!=null){
            kernelAlpha = spatialKernal.geta();
            this.addModel(spatialKernal);
        }
        this.transmissionRate = transmissionRate;
        hasLatentPeriods = outbreak.hasLatentPeriods();
        this.addModel(treeLikelihood);
        this.addVariable(transmissionRate);
        infectionTimes = treeLikelihood.getInfTimesMap();
        probFunct = new InnerIntegral();
        likelihoodKnown = false;
        hasGeography = spatialKernal!=null;
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if(model instanceof CaseToCaseTreeLikelihood){
            infectionTimes = treeLikelihood.getInfTimesMap();
        }
        likelihoodKnown = false;
    }

    protected void storeState() {
        storedLogLikelihood = logLikelihood;
        storedLikelihoodKnown = likelihoodKnown;
        storedInfectionTimes = new HashMap<Integer, Double>(infectionTimes);
    }

    protected void restoreState() {
        logLikelihood = storedLogLikelihood;
        likelihoodKnown = storedLikelihoodKnown;
        infectionTimes = storedInfectionTimes;
    }

    protected void acceptState() {
        // nothing to do
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        likelihoodKnown = false;
    }

    public Model getModel() {
        return this;
    }

    public double getLogLikelihood() {
        if(!likelihoodKnown){
            int N = outbreak.size();
            double lambda = transmissionRate.getParameterValue(0);
            double logUnnormalisedValue = Math.log(Math.pow(lambda, N-1)) + Math.log(aAlpha()) - lambda*bAlpha() - Math.log(N);
            double logNormalisationValue;
            if(hasGeography){
                // @todo don't calculate this if neither alpha nor the infection times have changed?
                logNormalisationValue = Math.log(probFunct.evaluateIntegral(0, kernelAlpha.getBounds().getUpperLimit(0)));
            } else {
                logNormalisationValue = Math.log(aAlpha()/Math.pow(bAlpha(), N));
            }
            logNormalisationValue += Math.log(MathUtils.factorial(N-1)) - Math.log(N);
            double treeLogL = treeLikelihood.getLogLikelihood();
            logLikelihood =  treeLogL + logUnnormalisedValue - logNormalisationValue;
            likelihoodKnown = true;
        }
        return logLikelihood;
    }


    public void makeDirty() {
        likelihoodKnown = false;
        treeLikelihood.makeDirty();
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

    // aAlpha for the current value of alpha

    private double aAlpha(){
        double product = 1;
        for(int i=1; i<outbreak.size(); i++){
            if(treeLikelihood.getInfector(i)!=null){
                product *= outbreak.getKernalValue(outbreak.getCase(i), treeLikelihood.getInfector(i), spatialKernel);
            }
        }
        return product;
    }

    // aAlpha for a given value of alpha

    private double aAlpha(double alpha){
        double product = 1;
        ArrayList<AbstractCase> copyOfCases = new ArrayList<AbstractCase>(outbreak.getCases());
        Collections.sort(copyOfCases, new CaseInfectionComparator());
        for(int i=1; i<outbreak.size(); i++){
            product *= outbreak.getKernalValue(copyOfCases.get(i), copyOfCases.get(i-1), spatialKernel, alpha);
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

        public double evaluateIntegral(double a, double b) {
            if(b!=Double.POSITIVE_INFINITY){
                return finiteIntegrator.integrate(this, a, b);
            } else {
                return infiniteIntegrator.integrate(this, a, b);
            }
        }

        public double evaluate(double argument) {
            return aAlpha(argument)/Math.pow(bAlpha(argument), outbreak.size());
        }

        public double getLowerBound() {
            return 0;
        }

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
                new ElementRule(SpatialKernel.class, "The spatial kernel", 0, 1),
                new ElementRule(TRANSMISSION_RATE, Parameter.class, "The transmission rate")
        };

    };



}


