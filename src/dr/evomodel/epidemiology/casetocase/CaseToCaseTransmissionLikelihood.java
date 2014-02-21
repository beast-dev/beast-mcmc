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

public class CaseToCaseTransmissionLikelihood extends AbstractModelLikelihood implements Loggable {

    private boolean DEBUG = false;

    private AbstractOutbreak outbreak;
    private CaseToCaseTreeLikelihood treeLikelihood;
    private SpatialKernel spatialKernel;
    private Parameter kernelAlpha;
    private Parameter transmissionRate;
    private boolean hasLatentPeriods;
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
    private final MultivariateIntegral integrator;
    private ArrayList<AbstractCase> sortedCases;
    private ArrayList<AbstractCase> storedSortedCases;
    private totalProbability totalProbability;

    private int integratorSteps;
    private int integratorBins;
    private double[] randomDoubleQueue;
    private double[] storedRandomDoubleQueue;

    public static final String CASE_TO_CASE_TRANSMISSION_LIKELIHOOD = "caseToCaseTransmissionLikelihood";

    public CaseToCaseTransmissionLikelihood(String name, AbstractOutbreak outbreak,
                                            CaseToCaseTreeLikelihood treeLikelihood, SpatialKernel spatialKernal,
                                            Parameter transmissionRate, int steps, int bins){
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
        this.integratorBins = bins;
        integrator = new QueuedMultivariateMonteCarloIntegral(integratorSteps, integratorBins);
        hasGeography = spatialKernal!=null;
        totalProbability = new totalProbability(hasGeography);
        hasLatentPeriods = treeLikelihood.hasLatentPeriods();
        sortCases();
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if(model instanceof CaseToCaseTreeLikelihood){
            // todo if you're going for maximum efficiency, some tree moves may not change the infection times, and most tree moves don't change MANY infection times
            treeProbKnown = false;
            if(!(object instanceof AbstractOutbreak)){
                transProbKnown = false;
                normalisationKnown = false;
                sortedCases = null;
                randomDoubleQueue = null;
            }
        } else if(model instanceof SpatialKernel){
            transProbKnown = false;
        }
        likelihoodKnown = false;
    }

    // no need to change the RNG queue unless the normalisation will have changed

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if(variable==transmissionRate){
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
        storedSortedCases = new ArrayList<AbstractCase>(sortedCases);
        storedRandomDoubleQueue = randomDoubleQueue;
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
        sortedCases = storedSortedCases;
        randomDoubleQueue = storedRandomDoubleQueue;
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
        integratorSteps = steps;
    }

    public void setIntegratorBins(int bins){
        this.integratorBins = bins;
    }

    private void testLoop(int gridSize, int tries, int maxSamples) throws IOException{
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

        double[][] values2 = new double[tries][maxSamples];

        for(int i=0; i<maxSamples; i++){
            setIntegratorSteps(i+1);
            for(int j=0; j<tries; j++){
                if(hasGeography){

                    randomDoubleQueue = new double[2*integratorSteps*integratorBins*integratorBins];
                    for(int k=0; k<2*integratorSteps*integratorBins*integratorBins; k++){
                        randomDoubleQueue[k] = MathUtils.nextDouble();
                    }


                    double expnormalisation = integrator.integrate(totalProbability,
                            new double[]{transmissionRate.getBounds().getLowerLimit(0),
                                    spatialKernel.geta().getBounds().getLowerLimit(0)},
                            new double[]{transmissionRate.getBounds().getUpperLimit(0),
                                    spatialKernel.geta().getBounds().getUpperLimit(0)});

                    values2[j][i] = Math.log(expnormalisation);

                } else {
                    if(randomDoubleQueue == null){
                        randomDoubleQueue = new double[integratorSteps];
                        for(int k=0; k<integratorSteps; k++){
                            randomDoubleQueue[k] = MathUtils.nextDouble();
                        }
                    }

                    double expnormalisation = integrator.integrate(totalProbability,
                            new double[]{transmissionRate.getBounds().getLowerLimit(0)},
                            new double[]{transmissionRate.getBounds().getUpperLimit(0)});

                    values2[j][i] = Math.log(expnormalisation);
                }

            }

        }

        for(int i=1; i<=maxSamples; i++){
            writer2.write(Integer.toString(i));
            if(i<maxSamples){
                writer2.write(",");
            }
        }
        writer2.newLine();

        for(int i=0; i<tries; i++){
            for(int j=0; j<maxSamples; j++){
                writer2.write(Double.toString(values2[i][j]));
                if(j<maxSamples){
                    writer2.write(",");
                }
            }
            writer2.newLine();
        }

        writer2.flush();
        writer2.close();

    }

    public double getLogLikelihood() {

        if(DEBUG){
            treeLikelihood.debugOutputTree("test.nex", false);
        }

        if(!likelihoodKnown){
            if(!treeProbKnown){
                treeLikelihood.prepareTimings();
            }
            if(!transProbKnown){
                try{
                    if(hasGeography){
                        transLogProb = Math.log(totalProbability.evaluate(new double[]
                                {transmissionRate.getParameterValue(0),spatialKernel.geta().getParameterValue(0)}));
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
                if(DEBUG){
                    try{
                        testLoop(100, 50, 5);
                    } catch(IOException e){
                        e.printStackTrace();
                    }
                }

                if(hasGeography){
                    if(randomDoubleQueue == null){
                        randomDoubleQueue = new double[2*(int)Math.pow(integratorBins,2)*integratorSteps];
                        for(int i=0; i<2*(int)Math.pow(integratorBins,2)*integratorSteps; i++){
                            randomDoubleQueue[i] = MathUtils.nextDouble();
                        }
                    }

                    double expNormalisation = integrator.integrate(totalProbability,
                            new double[]{transmissionRate.getBounds().getLowerLimit(0),
                                    spatialKernel.geta().getBounds().getLowerLimit(0)},
                            new double[]{transmissionRate.getBounds().getUpperLimit(0),
                                    spatialKernel.geta().getBounds().getUpperLimit(0)});

                    normalisation = Math.log(expNormalisation);

                } else {
                    if(randomDoubleQueue == null){
                        randomDoubleQueue = new double[(int)Math.pow(integratorBins,2)*integratorSteps];
                        for(int i=0; i<(int)Math.pow(integratorBins,2)*integratorSteps; i++){
                            randomDoubleQueue[i] = MathUtils.nextDouble();
                        }
                    }

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



    private class totalProbability implements MultivariateFunction {

        final boolean hasGeography;

        private totalProbability(boolean hasGeography){
            this.hasGeography = hasGeography;
        }

        // index 0 is lambda, index 1 if present is alpha

        public double evaluate(double[] argument) {
            if(sortedCases == null){
                sortCases();
            }
            double logProb = 0;

            for(AbstractCase aCase : outbreak.getCases()){
                logProb += caseLogProbability(aCase, argument);
            }

            return Math.exp(logProb);
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
        public static final String INTEGRATOR_BINS = "integratorBins";

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            CaseToCaseTreeLikelihood c2cTL = (CaseToCaseTreeLikelihood)
                    xo.getChild(CaseToCaseTreeLikelihood.class);
            SpatialKernel kernel = (SpatialKernel) xo.getChild(SpatialKernel.class);
            Parameter transmissionRate = (Parameter) xo.getElementFirstChild(TRANSMISSION_RATE);
            int steps = 2;

            if(xo.hasAttribute(INTEGRATOR_STEPS)){
                steps = Integer.parseInt((String)xo.getAttribute(INTEGRATOR_STEPS));
            }

            int bins = 4;

            if(xo.hasAttribute(INTEGRATOR_BINS)){
                bins = Integer.parseInt((String)xo.getAttribute(INTEGRATOR_BINS));
            }


            return new CaseToCaseTransmissionLikelihood(CASE_TO_CASE_TRANSMISSION_LIKELIHOOD, c2cTL.getOutbreak(),
                    c2cTL, kernel, transmissionRate, steps, bins);
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
                AttributeRule.newIntegerRule(INTEGRATOR_STEPS, true),
                AttributeRule.newIntegerRule(INTEGRATOR_BINS, true)
        };

    };

    // Not the most elegant solution, but you want two types of log out of this model, one for numerical parameters
    // (which Tracer can read) and one for the transmission tree (which it cannot). This is set up so that C2CTransL
    // is the numerical log and C2CTreeL the TT one.

    public LogColumn[] getColumns(){
        ArrayList<LogColumn> columns = new ArrayList<LogColumn>();
        for(int i=0; i<outbreak.size(); i++){
            final AbstractCase infected = outbreak.getCase(i);
            columns.add(new LogColumn.Abstract(infected.toString()+"_infection_date"){
                protected String getFormattedValue() {
                    return String.valueOf(treeLikelihood.getInfectionTime(infected));
                }
            });
            if(hasLatentPeriods){
                columns.add(new LogColumn.Abstract(infected.toString()+"_infectious_date"){
                    protected String getFormattedValue() {
                        return String.valueOf(treeLikelihood.getInfectiousTime(infected));
                    }
                });
                columns.add(new LogColumn.Abstract(infected.toString()+"_latent_period"){
                    protected String getFormattedValue() {
                        return String.valueOf(treeLikelihood.getLatentPeriod(infected));
                    }
                });
            }
            columns.add(new LogColumn.Abstract(infected.toString()+"_infectious_period"){
                protected String getFormattedValue() {
                    return String.valueOf(treeLikelihood.getInfectiousPeriod(infected));
                }
            });
            if(hasLatentPeriods){
                columns.add(new LogColumn.Abstract(infected.toString()+"_infected_period"){
                    protected String getFormattedValue() {
                        return String.valueOf(
                                treeLikelihood.getInfectiousPeriod(infected)+treeLikelihood.getLatentPeriod(infected));
                    }
                });
            }
        }
        columns.add(new LogColumn.Abstract("infectious_period.mean"){
            protected String getFormattedValue() {
                return String.valueOf(CaseToCaseTreeLikelihood
                        .getSummaryStatistics(treeLikelihood.getInfectiousPeriods(false))[0]);
            }
        });
        columns.add(new LogColumn.Abstract("infectious_period.median"){
            protected String getFormattedValue() {
                return String.valueOf(CaseToCaseTreeLikelihood
                        .getSummaryStatistics(treeLikelihood.getInfectiousPeriods(false))[1]);
            }
        });
        columns.add(new LogColumn.Abstract("infectious_period.var") {
            protected String getFormattedValue() {
                return String.valueOf(CaseToCaseTreeLikelihood
                        .getSummaryStatistics(treeLikelihood.getInfectiousPeriods(false))[2]);
            }
        });
        columns.add(new LogColumn.Abstract("infectious_period.stdev"){
            protected String getFormattedValue() {
                return String.valueOf(CaseToCaseTreeLikelihood
                        .getSummaryStatistics(treeLikelihood.getInfectiousPeriods(false))[3]);
            }
        });
        if(hasLatentPeriods){
            columns.add(new LogColumn.Abstract("latent_period.mean"){
                protected String getFormattedValue() {
                    return String.valueOf(CaseToCaseTreeLikelihood
                            .getSummaryStatistics(treeLikelihood.getLatentPeriods(false))[0]);
                }
            });
            columns.add(new LogColumn.Abstract("latent_period.median"){
                protected String getFormattedValue() {
                    return String.valueOf(CaseToCaseTreeLikelihood
                            .getSummaryStatistics(treeLikelihood.getLatentPeriods(false))[1]);
                }
            });
            columns.add(new LogColumn.Abstract("latent_period.var") {
                protected String getFormattedValue() {
                    return String.valueOf(CaseToCaseTreeLikelihood
                            .getSummaryStatistics(treeLikelihood.getLatentPeriods(false))[2]);
                }
            });
            columns.add(new LogColumn.Abstract("latent_period.stdev"){
                protected String getFormattedValue() {
                    return String.valueOf(CaseToCaseTreeLikelihood
                            .getSummaryStatistics(treeLikelihood.getLatentPeriods(false))[3]);
                }
            });
            columns.add(new LogColumn.Abstract("infected_period.mean"){
                protected String getFormattedValue() {
                    return String.valueOf(CaseToCaseTreeLikelihood
                            .getSummaryStatistics(treeLikelihood.getInfectedPeriods(false))[0]);
                }
            });
            columns.add(new LogColumn.Abstract("infected_period.median"){
                protected String getFormattedValue() {
                    return String.valueOf(CaseToCaseTreeLikelihood
                            .getSummaryStatistics(treeLikelihood.getInfectedPeriods(false))[1]);
                }
            });
            columns.add(new LogColumn.Abstract("infected_period.var") {
                protected String getFormattedValue() {
                    return String.valueOf(CaseToCaseTreeLikelihood
                            .getSummaryStatistics(treeLikelihood.getInfectedPeriods(false))[2]);
                }
            });
            columns.add(new LogColumn.Abstract("infected_period.stdev"){
                protected String getFormattedValue() {
                    return String.valueOf(CaseToCaseTreeLikelihood
                            .getSummaryStatistics(treeLikelihood.getInfectedPeriods(false))[3]);
                }
            });
        }

        return columns.toArray(new LogColumn[columns.size()]);
    }

    /*
    If an individual likelihood calculation involves MC integration then there is an obvious problem with full
    evaluation errors. To deal with this we need to queue up RNG draws beforehand and save them, clearing the queue
    on model change events that involve recalculating the normalisation but _not_ on makeDirty().
    */

    private class QueuedMultivariateMonteCarloIntegral extends MultivariateMonteCarloIntegral{

        private double[][] corners;

        public QueuedMultivariateMonteCarloIntegral(int sampleSize) {
            super(sampleSize);
        }

        public QueuedMultivariateMonteCarloIntegral(int sampleSize, int bins) {
            super(sampleSize, bins);
        }

        public double integrate(MultivariateFunction f, double[] mins, double[] maxes) {

            int dim = f.getNumArguments();
            int totalBins = (int)Math.pow(integratorBins, dim);
            double[] steps = new double[dim];
            double totalArea=1;


            for(int i=0; i<dim; i++){
                totalArea *= (maxes[i]-mins[i]);
                steps[i] = (maxes[i]-mins[i])/integratorBins;
            }

            if(corners==null){

                corners = new double[totalBins][dim];
                int[] currentGridPosition = new int[dim];


                for(int index=0; index<totalBins; index++){
                    double[] currentCorner = new double[dim];
                    for(int i=0; i<dim; i++){
                        currentCorner[i] = mins[i] + currentGridPosition[i]*steps[i];
                    }

                    corners[index]=currentCorner;

                    int dimToCheck = 0;
                    while(dimToCheck<dim){
                        if(currentGridPosition[dimToCheck]+1<getBins()){
                            currentGridPosition[dimToCheck]++;
                            break;
                        } else {
                            currentGridPosition[dimToCheck] = 0;
                        }
                        dimToCheck++;
                    }
                }
            }

            double integral = 0.0;

            int queuePosition = 0;

            for(int i=0; i<totalBins; i++){

                for (int j=1; j <= integratorSteps; j++) {

                    double[] sample = new double[dim];
                    for(int k=0; k<sample.length; k++){
                        sample[k] = corners[i][k] + randomDoubleQueue[queuePosition]*(steps[k]);
                        queuePosition++;
                    }

                    integral += f.evaluate(sample);
                }

            }
            integral *= totalArea/((double)integratorSteps*totalBins);
            return integral;
        }

    }

}
