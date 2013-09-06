package dr.evomodel.epidemiology.casetocase;

import dr.inference.model.*;
import dr.math.IntegrableUnivariateFunction;
import dr.math.RiemannApproximation;
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
    private HashMap<Integer, Double> infectionTimes;
    private HashMap<Integer, Double> storedInfectionTimes;
    private IntegrableUnivariateFunction probFunct;
    private boolean likelihoodKnown;
    private boolean storedLikelihoodKnown;
    private boolean transProbKnown;
    private boolean storedTransProbKnown;
    private boolean normalisationKnown;
    private boolean storedNormalisationKnown;
    private boolean treeProbKnown;
    private boolean storedTreeProbKnown;
    private double logLikelihood;
    private double storedLogLikelihood;
    private double transLogProb;
    private double storedTransLogProb;
    private double normalisation;
    private double storedNormalisation;
    private double treeLogProb;
    private double storedTreeLogProb;
    private boolean hasGeography;
    private boolean integrateToInfinity;
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
            integrateToInfinity =  kernelAlpha.getBounds().getUpperLimit(0)==Double.POSITIVE_INFINITY;
        }
        this.transmissionRate = transmissionRate;
        hasLatentPeriods = outbreak.hasLatentPeriods();
        this.addModel(treeLikelihood);
        this.addVariable(transmissionRate);
        infectionTimes = treeLikelihood.getInfTimesMap();
        if(spatialKernal!=null && integrateToInfinity){
            probFunct = new InnerIntegralTransformed(new InnerIntegral());
        } else {
            probFunct = new InnerIntegral();
        }
        likelihoodKnown = false;
        hasGeography = spatialKernal!=null;
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if(model instanceof CaseToCaseTreeLikelihood){
            // @todo if you're going for maximum efficiency, some tree moves may not change the infection times, and most tree moves don't change MANY infection times
            treeProbKnown = false;
            transProbKnown = false;
            normalisationKnown = false;
        } else if(model instanceof SpatialKernel){
            transProbKnown = false;
        }
        likelihoodKnown = false;
    }

    protected void storeState() {
        storedLogLikelihood = logLikelihood;
        storedLikelihoodKnown = likelihoodKnown;
        storedNormalisation = normalisation;
        storedNormalisationKnown = normalisationKnown;
        storedTransLogProb = transLogProb;
        storedTransProbKnown = transProbKnown;
        storedTreeLogProb = treeLogProb;
        storedTreeProbKnown = treeProbKnown;
        storedInfectionTimes = new HashMap<Integer, Double>(infectionTimes);
    }

    protected void restoreState() {
        logLikelihood = storedLogLikelihood;
        likelihoodKnown = storedLikelihoodKnown;
        transLogProb = storedTransLogProb;
        transProbKnown = storedTransProbKnown;
        treeLogProb = storedTreeLogProb;
        treeProbKnown = storedTreeProbKnown;
        normalisation = storedNormalisation;
        normalisationKnown = storedNormalisationKnown;
        infectionTimes = storedInfectionTimes;
    }

    protected void acceptState() {
        // nothing to do
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if(variable==transmissionRate){
            transProbKnown = false;
        }
        likelihoodKnown = false;
    }

    public Model getModel() {
        return this;
    }

    public double getLogLikelihood() {
        if(!likelihoodKnown){
            int N = outbreak.size();
            double lambda = transmissionRate.getParameterValue(0);
            if(!treeProbKnown){
                infectionTimes = treeLikelihood.getInfTimesMap();
            }
            if(!transProbKnown){
                transLogProb = Math.log(Math.pow(lambda, N-1)) + Math.log(aAlpha()) - lambda*bAlpha() - Math.log(N);
            }
            if(!normalisationKnown){
                if(hasGeography){
                    if(integrateToInfinity){
                        normalisation = Math.log(probFunct.evaluateIntegral(0, 1));
                    } else {
                        normalisation = Math.log(probFunct.evaluateIntegral(0, kernelAlpha.getBounds().getUpperLimit(0)));
                    }
                } else {
                    normalisation = Math.log(aAlpha()/Math.pow(bAlpha(), N));
                }
                normalisation += Math.log(MathUtils.factorial(N-1)) - Math.log(N);

            }
            if(!treeProbKnown){
                treeLogProb = treeLikelihood.getLogLikelihood();
            }
            logLikelihood =  treeLogProb + transLogProb - normalisation;
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

        RiemannApproximation integrator;

        private InnerIntegral(){
            integrator = new RiemannApproximation(50);
        }

        public double evaluateIntegral(double a, double b) {
            return integrator.integrate(this, a, b);
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

    // integral of the former from 0 to +Infinity

    private class InnerIntegralTransformed implements IntegrableUnivariateFunction {

        RiemannApproximation integrator;
        InnerIntegral originalFunction;

        private InnerIntegralTransformed(InnerIntegral inner){
            this.originalFunction = inner;
            integrator = new RiemannApproximation(50);
        }

        public double evaluateIntegral(double a, double b) {
            return integrator.integrate(this, a, b);
        }

        public double evaluate(double argument) {
            return originalFunction.evaluate(argument/(1-Math.pow(argument,2)))/Math.pow(1-argument,2);
        }

        public double getLowerBound() {
            return 0;
        }

        public double getUpperBound() {
            return 1;
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


