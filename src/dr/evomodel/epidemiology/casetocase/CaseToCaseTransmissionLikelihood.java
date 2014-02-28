package dr.evomodel.epidemiology.casetocase;

import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.model.*;
import dr.math.*;
import dr.xml.*;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

/**
 * A likelihood function for transmission between identified epidemiological outbreak
 *
 * Timescale must be in days. Python scripts to write XML for it and analyse the posterior set of networks exist;
 * contact MH. @todo make timescale not just in days
 *
 * Latent periods are not implemented currently
 *
 * @author Matthew Hall
 * @version $Id: $
 */

public class CaseToCaseTransmissionLikelihood extends AbstractModelLikelihood implements Loggable {

    private static final boolean DEBUG = true;
    private static final double INTEGRAL_APPROXIMATION_TOLERANCE = 1;
    private static final int RETEST_INTEGRAL = 1000;
    private static final int FULL_EVALUATION = 10000;
    private static final int MAX_STEPS=10;
    private static final int MIN_STEPS=5;
    private int mcmcState = -1;

    private AbstractOutbreak outbreak;
    private CaseToCaseTreeLikelihood treeLikelihood;
    private SpatialKernel spatialKernel;
    private Parameter kernelAlpha;
    private Parameter transmissionRate;
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
    private final boolean hasGeography;
    private final BespokeNumericalIntegrator integrator;
    private ArrayList<AbstractCase> sortedCases;
    private ArrayList<AbstractCase> storedSortedCases;
    private TotalProbability totalProbability;

    private int integratorSteps;
    private int storedIntegratorSteps;

    public static final String CASE_TO_CASE_TRANSMISSION_LIKELIHOOD = "caseToCaseTransmissionLikelihood";

    public CaseToCaseTransmissionLikelihood(String name, AbstractOutbreak outbreak,
                                            CaseToCaseTreeLikelihood treeLikelihood, SpatialKernel spatialKernal,
                                            Parameter transmissionRate, int steps){
        super(name);
        this.outbreak = outbreak;
        this.treeLikelihood = treeLikelihood;
        this.spatialKernel = spatialKernal;
        if(spatialKernal!=null){
            kernelAlpha = spatialKernal.geta();
            if(kernelAlpha.getBounds().getUpperLimit(0)==Double.POSITIVE_INFINITY){
                throw new RuntimeException("Infinite upper limits for kernel parameters not implemeted yet");
            }
            this.addModel(spatialKernal);
        }
        this.transmissionRate = transmissionRate;
        if(transmissionRate.getBounds().getUpperLimit(0)==Double.POSITIVE_INFINITY){
            throw new RuntimeException("Infinite upper limits for transmission rate not implemented yet");
        }
        this.addModel(treeLikelihood);
        this.addVariable(transmissionRate);
        likelihoodKnown = false;
        this.integratorSteps = steps;
        integrator = new BespokeNumericalIntegrator(integratorSteps);
        hasGeography = spatialKernal!=null;
        totalProbability = new TotalProbability(hasGeography);
        sortCases();
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if(model instanceof CaseToCaseTreeLikelihood){

            treeProbKnown = false;
            if(!(object instanceof AbstractOutbreak)){
                transProbKnown = false;
                normalisationKnown = false;
                sortedCases = null;
            }
        } else if(model instanceof SpatialKernel){
            transProbKnown = false;

            // because of the way the numerical integrator works currently, changing the parameter does change
            // the normalisation

            normalisationKnown = false;
        }
        likelihoodKnown = false;
    }

    // no need to change the RNG queue unless the normalisation will have changed

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if(variable==transmissionRate){
            transProbKnown = false;

            // because of the way the numerical integrator works currently, changing the parameter does change
            // the normalisation

            normalisationKnown = false;
        }
        likelihoodKnown = false;
    }

    protected void storeState() {
        // I apologise for this horrible hack. The number of times storeState() is has been called, minus one, is
        // the current state number of the chain, and the dynamic monitoring of the numerical integrator performance
        // needs to know the current state number
        mcmcState++;


        storedLogLikelihood = logLikelihood;
        storedLikelihoodKnown = likelihoodKnown;
        storedNormalisation = normalisation;
        storedNormalisationKnown = normalisationKnown;
        storedTransLogProb = transLogProb;
        storedTransProbKnown = transProbKnown;
        storedTreeLogProb = treeLogProb;
        storedTreeProbKnown = treeProbKnown;
        storedIntegratorSteps = integrator.steps;
        storedSortedCases = new ArrayList<AbstractCase>(sortedCases);
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
        integrator.steps = storedIntegratorSteps;
        sortedCases = storedSortedCases;
    }

    protected void acceptState() {
        // nothing to do
    }

    public Model getModel() {
        return this;
    }

    public CaseToCaseTreeLikelihood getTreeLikelihood(){
        return treeLikelihood;
    }

    public void setIntegratorSteps(int steps){
        integrator.steps = steps;
    }

    private void testLoop(int gridSize, int maxSamples) throws IOException{
        BufferedWriter writer = new BufferedWriter(new FileWriter("forgraph.csv"));

        double[] steps = new double[2];
        steps[0] = (transmissionRate.getBounds().getUpperLimit(0)-  transmissionRate.getBounds().getLowerLimit(0))
                /gridSize;
        steps[1] = (spatialKernel.geta().getBounds().getUpperLimit(0)-spatialKernel.geta().getBounds().getLowerLimit(0))
                /gridSize;

        writer.write(",");
        for(int i=0; i<=gridSize; i++){
            writer.write(Double.toString(i*steps[1]));
            if(i<gridSize){
                writer.write(",");
            }
        }
        writer.newLine();

        for(int i=0; i<=gridSize; i++){
            writer.write(Double.toString(i*steps[0])+",");
            for(int j=0; j<=gridSize; j++){
                writer.write(Double.toString(totalProbability.evaluate(new double[]
                        {i*steps[0], j*steps[1]})));
                if(j<gridSize){
                    writer.write(",");
                }
            }
            writer.newLine();
        }

        writer.flush();
        writer.close();


        BufferedWriter writer2 = new BufferedWriter(new FileWriter("MCvartest2.csv"));

        double[] values2 = new double[maxSamples];

        for(int i=0; i<maxSamples; i++){
            setIntegratorSteps(i+1);

            if(hasGeography){

                double expnormalisation = integrator.integrate(totalProbability,
                        new double[]{transmissionRate.getBounds().getLowerLimit(0),
                                spatialKernel.geta().getBounds().getLowerLimit(0)},
                        new double[]{transmissionRate.getBounds().getUpperLimit(0),
                                spatialKernel.geta().getBounds().getUpperLimit(0)});

                values2[i] = Math.log(expnormalisation);

            } else {

                double expnormalisation = integrator.integrate(totalProbability,
                        new double[]{transmissionRate.getBounds().getLowerLimit(0)},
                        new double[]{transmissionRate.getBounds().getUpperLimit(0)});

                values2[i] = Math.log(expnormalisation);
            }



        }

        for(int i=1; i<=maxSamples; i++){
            writer2.write(Integer.toString(i));
            if(i<maxSamples){
                writer2.write(",");
            }
        }
        writer2.newLine();


        for(int j=0; j<maxSamples; j++){
            writer2.write(Double.toString(values2[j]));
            if(j<maxSamples){
                writer2.write(",");
            }
        }
        writer2.newLine();


        writer2.flush();
        writer2.close();

    }

    public double getLogLikelihood() {

//        if(DEBUG){
//            treeLikelihood.debugOutputTree("test.nex", false);
//        }

        if(!likelihoodKnown){
            if(!treeProbKnown){
                treeLikelihood.prepareTimings();
            }
            if(!transProbKnown){
                try{
                    if(hasGeography){
                        double[] point = {transmissionRate.getParameterValue(0),
                                spatialKernel.geta().getParameterValue(0)};

                        transLogProb = Math.log(totalProbability.evaluate(point));
                    } else {
                        transLogProb = Math.log(totalProbability.evaluate(new double[]
                                {transmissionRate.getParameterValue(0)}));
                    }
                } catch(BadPartitionException e){
                    transLogProb = Double.NEGATIVE_INFINITY;
                    logLikelihood = Double.NEGATIVE_INFINITY;
                    transProbKnown = true;
                    likelihoodKnown = true;
                    return logLikelihood;
                }
                transProbKnown = true;
            }
            if(!normalisationKnown){
//                if(DEBUG){
//                    try{
//                        testLoop(100, 20);
//                    } catch(IOException e){
//                        e.printStackTrace();
//                    }
//                }
                if(hasGeography){

                    double[] point;
                    point = new double[]{transmissionRate.getParameterValue(0),
                            spatialKernel.geta().getParameterValue(0)};

                    integrator.centre = point;

                    normalisation = integrator.logIntegrate(totalProbability,
                            new double[]{transmissionRate.getBounds().getLowerLimit(0),
                                    spatialKernel.geta().getBounds().getLowerLimit(0)},
                            new double[]{transmissionRate.getBounds().getUpperLimit(0),
                                    spatialKernel.geta().getBounds().getUpperLimit(0)});

                    while(normalisation == Double.NEGATIVE_INFINITY){
                        if(integrator.steps<MAX_STEPS){
                            setIntegratorSteps(integrator.steps+1);
                            if(DEBUG){
                                System.out.println("Increasing the number of integrator steps to "+integrator.steps
                                        +" due to underflow");
                            }
                            normalisation = integrator.logIntegrate(totalProbability,
                                    new double[]{transmissionRate.getBounds().getLowerLimit(0),
                                            spatialKernel.geta().getBounds().getLowerLimit(0)},
                                    new double[]{transmissionRate.getBounds().getUpperLimit(0),
                                            spatialKernel.geta().getBounds().getUpperLimit(0)});
                        } else {
                            throw new RuntimeException("Infinite normalisation but can't increase number of steps " +
                                    "further");
                        }
                    }

                    // periodically test that the numerical integrator is working OK
                    // we're assuming that increasing the sample size always increases accuracy.

                    if(mcmcState<=FULL_EVALUATION || mcmcState % RETEST_INTEGRAL == 0){

                        int baseSteps = integrator.steps;

                        boolean changedSteps = false;

                        int currentSteps = integrator.steps;

                        double testSlower = normalisation;

                        boolean ok = false;

                        while(!ok && currentSteps<=MAX_STEPS){
                            double old = testSlower;

                            currentSteps = integrator.steps;

                            setIntegratorSteps(currentSteps+1);

                            testSlower = integrator.logIntegrate(totalProbability,
                                    new double[]{transmissionRate.getBounds().getLowerLimit(0),
                                            spatialKernel.geta().getBounds().getLowerLimit(0)},
                                    new double[]{transmissionRate.getBounds().getUpperLimit(0),
                                            spatialKernel.geta().getBounds().getUpperLimit(0)});

                            if(Math.abs(testSlower-old)<INTEGRAL_APPROXIMATION_TOLERANCE){
                                ok = true;
                                setIntegratorSteps(currentSteps);

                                if(currentSteps!=baseSteps){
                                    changedSteps = true;
                                    if(DEBUG){
                                        System.out.println("Increasing the number of integrator steps to "
                                                +integrator.steps+" due to unreliability");
                                    }
                                }

                                normalisation = old;


                            } else if(currentSteps==MAX_STEPS){
                                ok = true;

                                setIntegratorSteps(currentSteps);

                                if(currentSteps!=baseSteps){
                                    changedSteps = true;
                                    if(DEBUG){
                                        System.out.println("Tried to increase the number of integrator steps to "
                                                +(integrator.steps+1)+" due to unreliability, but hit maximum " +
                                                "(difference="+Math.abs(testSlower-old)+"); stayed at "+MAX_STEPS);
                                    }
                                }

                                normalisation = old;


                            }

                        }
                        if(!changedSteps) {

                            // don't do this in the full evaluation phase, just to be safe

                            if(mcmcState>FULL_EVALUATION){

                                currentSteps = integrator.steps;

                                double testFaster = normalisation;

                                ok = true;

                                while(ok && currentSteps>=MIN_STEPS){

                                    double old = testFaster;

                                    currentSteps = integrator.steps;

                                    setIntegratorSteps(currentSteps-1);

                                    testFaster = integrator.logIntegrate(totalProbability,
                                            new double[]{transmissionRate.getBounds().getLowerLimit(0),
                                                    spatialKernel.geta().getBounds().getLowerLimit(0)},
                                            new double[]{transmissionRate.getBounds().getUpperLimit(0),
                                                    spatialKernel.geta().getBounds().getUpperLimit(0)});


                                    if(Math.abs(testFaster-old)>INTEGRAL_APPROXIMATION_TOLERANCE){
                                        ok = false;
                                        setIntegratorSteps(currentSteps);

                                        normalisation = old;

                                        if(currentSteps!=baseSteps){
                                            changedSteps = true;
                                            if(DEBUG){
                                                System.out.println("Decreasing the number of integrator steps to "
                                                        +integrator.steps);
                                            }
                                        }


                                    } else if(currentSteps==MIN_STEPS){
                                        ok = false;
                                        setIntegratorSteps(currentSteps);

                                        normalisation = old;

                                        if(currentSteps!=baseSteps){
                                            changedSteps = true;
                                            if(DEBUG){
                                                System.out.println("Tried to decrease the number of integrator steps " +
                                                        "to "+(integrator.steps-1)+", but hit minimum (distance="
                                                        +Math.abs(testFaster-old)+"); stayed at "+MIN_STEPS);
                                            }
                                        }



                                    }
                                }

                            }
                            if(!changedSteps){
                                setIntegratorSteps(baseSteps);
                            }
                        }
                    }


                } else {

                    //todo this needs updating

                    double expNormalisation = integrator.integrate(totalProbability,
                            new double[]{transmissionRate.getBounds().getLowerLimit(0)},
                            new double[]{transmissionRate.getBounds().getUpperLimit(0)});

                    normalisation = Math.log(expNormalisation);
                }

                normalisationKnown = true;
            }
            if(!treeProbKnown){
                treeLogProb = treeLikelihood.getLogLikelihood();
                treeProbKnown = true;
            }
            logLikelihood =  treeLogProb + transLogProb - normalisation;
            likelihoodKnown = true;
        }

        return logLikelihood;
    }


    // Gibbs operator needs this

    public double calculateTempLogLikelihood(AbstractCase[] map){

        // todo probably this should tell PartitionedTreeModel what needs recalculating

        BranchMapModel branchMap = treeLikelihood.getBranchMap();

        AbstractCase[] trueMap = branchMap.getArrayCopy();
        branchMap.setAll(map, false);
        double out = getLogLikelihood();
        branchMap.setAll(trueMap, false);

        return out;
    }


    public void makeDirty() {
        likelihoodKnown = false;
        transProbKnown = false;
        treeProbKnown = false;
        normalisationKnown = false;
        sortedCases = null;
        treeLikelihood.makeDirty();
    }

    private class CaseInfectionComparator implements Comparator<AbstractCase> {
        public int compare(AbstractCase abstractCase, AbstractCase abstractCase2) {
            return Double.compare(treeLikelihood.getInfectionTime(abstractCase),
                    treeLikelihood.getInfectionTime(abstractCase2));
        }
    }


    private void sortCases(){
        sortedCases = new ArrayList<AbstractCase>(outbreak.getCases());
        Collections.sort(sortedCases, new CaseInfectionComparator());
    }

    private double caseLogProbability(AbstractCase infectee, double[] argument){
        double logP = 0;

        if(sortedCases == null){
            sortCases();
        }
        if(treeLikelihood.getInfector(infectee)==null){
            return -Math.log(outbreak.size());
        }
        AbstractCase infector = treeLikelihood.getInfector(infectee);

        // probability that the infector infected this case at this time

        double r = argument[0];

        if(argument.length>1){
            double kernelValue = outbreak.getKernelValue(infectee, infector, spatialKernel, argument[1]);
            r *= kernelValue;
        }

        double infecteeInfected = treeLikelihood.getInfectionTime(infectee);
        double infectorInfectious = treeLikelihood.getInfectiousTime(infector);
        double infecteeNoninfectious = infectee.getCullTime();
        double infectorNoninfectious = infector.getCullTime();

        // need to differentiate between an actual zero probability and a rounding error. This is a zero probability.

        if(infecteeInfected < infectorInfectious || infecteeInfected > infectorNoninfectious){
            throw new BadPartitionException("Illegal partition given known timings");
        }

        logP += Math.log(r);

        // probability that nothing infected this case sooner

        for(AbstractCase nonInfector : sortedCases){
            if(nonInfector == infectee){
                // no need to consider any infections happening after this one
                break;
            }

            r = argument[0];
            if(argument.length>1){
                double kernelValue = outbreak.getKernelValue(infectee, nonInfector, spatialKernel, argument[1]);
                r *= kernelValue;
            }
            double nonInfectorInfectious = treeLikelihood.getInfectiousTime(nonInfector);
            if(nonInfectorInfectious <= infecteeInfected){
                logP += -r*(infecteeInfected - nonInfectorInfectious);
            }
        }

        // probability that _something_ infected this case before it became noninfectious

        double logProduct = 0;

        for(AbstractCase nonInfector : sortedCases){
            if(nonInfector != infectee){
                double nonInfectorInfected = treeLikelihood.getInfectionTime(nonInfector);

                if(nonInfectorInfected > infecteeNoninfectious){
                    // no need to consider any infections happening later
                    break;
                }

                r = argument[0];
                if(argument.length>1){
                    double kernelValue = outbreak.getKernelValue(infectee, nonInfector, spatialKernel, argument[1]);
                    r *= kernelValue;
                }
                double nonInfectorInfectious = treeLikelihood.getInfectiousTime(nonInfector);
                if(nonInfectorInfectious <= infecteeNoninfectious){
                    logProduct += -r*(infecteeNoninfectious - nonInfectorInfectious);
                }
            }
        }

        double product = Math.exp(logProduct);

        logP -= Math.log1p(-product);

        return logP;
    }

    private class TotalProbability implements MultivariateFunction {

        final boolean hasGeography;

        private TotalProbability(boolean hasGeography){
            this.hasGeography = hasGeography;
        }

        // index 0 is lambda, index 1 if present is alpha

        public double evaluate(double[] argument) {
            return Math.exp(logEvaluate(argument));
        }

        public double logEvaluate(double[] argument) {
            if(sortedCases == null){
                sortCases();
            }
            double logProb = 0;

            for(AbstractCase aCase : outbreak.getCases()){
                logProb += caseLogProbability(aCase, argument);
            }

            return logProb;
        }

        public int getNumArguments() {
            return 2;
        }

        public double getLowerBound(int n) {
            return 0;
        }

        public double getUpperBound(int n) {
            return n==0 ? kernelAlpha.getBounds().getUpperLimit(0) : transmissionRate.getBounds().getUpperLimit(0);
        }
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public static final String TRANSMISSION_RATE = "transmissionRate";
        public static final String INTEGRATOR_STEPS = "integratorSteps";

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            CaseToCaseTreeLikelihood c2cTL = (CaseToCaseTreeLikelihood)
                    xo.getChild(CaseToCaseTreeLikelihood.class);
            SpatialKernel kernel = (SpatialKernel) xo.getChild(SpatialKernel.class);
            Parameter transmissionRate = (Parameter) xo.getElementFirstChild(TRANSMISSION_RATE);
            int steps = 2;

            if(xo.hasAttribute(INTEGRATOR_STEPS)){
                steps = Integer.parseInt((String)xo.getAttribute(INTEGRATOR_STEPS));
            }

            return new CaseToCaseTransmissionLikelihood(CASE_TO_CASE_TRANSMISSION_LIKELIHOOD, c2cTL.getOutbreak(),
                    c2cTL, kernel, transmissionRate, steps);
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
                new ElementRule(TRANSMISSION_RATE, Parameter.class, "The transmission rate"),
                AttributeRule.newIntegerRule(INTEGRATOR_STEPS, true)
        };

    };

    // Not the most elegant solution, but you want two types of log out of this model, one for numerical parameters
    // (which Tracer can read) and one for the transmission tree (which it cannot). This is set up so that C2CTransL
    // is the numerical log and C2CTreeL the TT one.

    public LogColumn[] getColumns(){
        return treeLikelihood.passColumns();
    }

    /*
    Awkwardly, we need a 2D numerical integrator to normalise the transmission probability. We have already sampled
    one nonzero point from this distribution; let's assume that's reasonably close to the mode in each dimension.

    Naive or stratified Monte Carlo integration has proved to be problematic. Importance sampling guided by the
    already sampled point might work.
    */

    private class BespokeNumericalIntegrator implements MultivariateIntegral{

        private int steps;
        private double[] centre;

        public BespokeNumericalIntegrator(int steps){
            this.steps = steps;
        }


        public double logIntegrate(MultivariateFunction f, double[] mins, double[] maxes) {

            double totalArea=1;
//            double[] stepSizes = new double[2];

            for(int i=0; i<2; i++){
                totalArea *= (maxes[i]-mins[i]);
//                stepSizes[i] = (maxes[i]-mins[i])/steps;
            }


            double[][] stepSizes = new double[2][2];

            for(int i=0; i<2; i++){
                stepSizes[i][0] = (centre[i] - mins[i])/Math.floor(steps/2);
                stepSizes[i][1] = (maxes[i] - centre[i])/Math.floor(steps/2);
            }

            double[] xSteps = new double[steps];
            double[] ySteps = new double[steps];

            for(int i=0; i<steps; i++){
                if(i<=(steps-1)/2){
                    xSteps[i] = mins[0] + (i+0.5)*stepSizes[0][0];
                    ySteps[i] = mins[1] + (i+0.5)*stepSizes[1][0];
                } else if(steps % 2!=0 && i==(steps-1)/2) {
                    xSteps[i] = centre[0];
                    ySteps[i] = centre[1];
                } else {
                    xSteps[i] = centre[0] + (i-Math.ceil(steps/2)+0.5)*stepSizes[0][1];
                    ySteps[i] = centre[1] + (i-Math.ceil(steps/2)+0.5)*stepSizes[1][1];
                }
            }

            double integral = 0.0;

            for(int i=0; i<steps; i++){
                for(int j=0; j<steps; j++){

                    double[] point = {xSteps[i], ySteps[j]};

                    integral += f.evaluate(point);

                }
            }

//            double integral = 0;
//
//            for(int i=0; i<steps; i++){
//                for(int j=0; j<steps; j++){
//
//                    double[] point = {mins[0] + (i+0.5)*stepSizes[0], mins[1] + (j+0.5)*stepSizes[1]};
//
//                    integral += f.evaluate(point);
//
//                }
//            }

            integral *= totalArea/(Math.pow(steps,2));

            if(integral>0){
                return Math.log(integral);
            } else {

                double logIntegral = Double.NEGATIVE_INFINITY;

                for(int i=0; i<steps; i++){
                    for(int j=0; j<steps; j++){

                        double[] point = {xSteps[i], ySteps[j]};

                        logIntegral = LogTricks.logSum(logIntegral, ((TotalProbability)f).logEvaluate(point));

                    }
                }


                logIntegral += Math.log(totalArea) - 2*Math.log(steps);

                return logIntegral;
//
//
//                for(int i=0; i<steps; i++){
//                    for(int j=0; j<steps; j++){
//
//                        double[] point = {mins[0] + (i+0.5)*stepSizes[0], mins[1] + (j+0.5)*stepSizes[1]};
//
//                        logIntegral = LogTricks.logSum(logIntegral, ((TotalProbability)f).logEvaluate(point));
//
//                    }
//                }
//
//
//                logIntegral += Math.log(totalArea) - 2*Math.log(steps);
//
//                return logIntegral;
            }
        }

        public double integrate(MultivariateFunction f, double[] mins, double[] maxes){
            return Math.exp(logIntegrate(f, mins, maxes));
        }
    }

}
