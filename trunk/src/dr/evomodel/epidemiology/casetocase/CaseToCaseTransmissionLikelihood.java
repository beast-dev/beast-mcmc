package dr.evomodel.epidemiology.casetocase;

import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.model.*;
import dr.math.*;
import dr.xml.*;
import org.apache.commons.math.util.MathUtils;
import org.netlib.util.doubleW;

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

    private final static boolean DEBUG = false;

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
    private totalLogProbability totalLogProbability;

    // the last time that the case corresponding to the first index could have infected the case corresponding to
    // the second

    private double[][] lastTimesToInfect;
    private double[][] storedLastTimesToInfect;

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
            this.addModel(spatialKernal);
        }
        this.transmissionRate = transmissionRate;
        this.addModel(treeLikelihood);
        this.addVariable(transmissionRate);
        likelihoodKnown = false;
        integrator = new MultivariateMonteCarloIntegral(steps);
        hasGeography = spatialKernal!=null;
        totalLogProbability = new totalLogProbability(hasGeography);
        hasLatentPeriods = treeLikelihood.hasLatentPeriods();
        // todo use BEAST classes for this, which means turning the infectious and latent period objects into Parameters
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
                lastTimesToInfect = null;
            }
        } else if(model instanceof SpatialKernel){
            transProbKnown = false;
        }
        likelihoodKnown = false;
    }

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

        storedLastTimesToInfect = new double[outbreak.size()][];
        for(int i = 0; i < outbreak.size(); i++){
            double[] aMatrix = lastTimesToInfect[i];
            int aLength = aMatrix.length;
            storedLastTimesToInfect[i] = Arrays.copyOf(aMatrix, aLength);
        }

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
        lastTimesToInfect = storedLastTimesToInfect;
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


    public double getLogLikelihood() {
        if(!likelihoodKnown){
            int N = outbreak.size();
            double lambda = transmissionRate.getParameterValue(0);
            if(!treeProbKnown){
                treeLikelihood.prepareTimings();
            }
            if(!transProbKnown){
                if(hasGeography){
                    transLogProb = totalLogProbability.evaluate(new double[]{transmissionRate.getParameterValue(0),
                            spatialKernel.geta().getParameterValue(0)});
                } else {
                    transLogProb = totalLogProbability.evaluate(new double[]{transmissionRate.getParameterValue(0)});
                }
                transProbKnown = true;
            }
            if(!normalisationKnown){
                if(hasGeography){
                    normalisation = integrator.integrate(totalLogProbability,
                            new double[]{transmissionRate.getBounds().getLowerLimit(0),
                                    spatialKernel.geta().getBounds().getLowerLimit(0)},
                            new double[]{transmissionRate.getBounds().getUpperLimit(0),
                                    spatialKernel.geta().getBounds().getUpperLimit(0)});

                } else {
                    normalisation = integrator.integrate(totalLogProbability,
                            new double[]{transmissionRate.getBounds().getLowerLimit(0)},
                            new double[]{spatialKernel.geta().getBounds().getLowerLimit(0)});
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
        lastTimesToInfect = null;
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

        if(infecteeInfected > infectorInfectious){
            return Double.NEGATIVE_INFINITY;
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
            double nonInfectorInfected = treeLikelihood.getInfectionTime(nonInfector);

            if(nonInfectorInfected > infecteeNoninfectious){
                // no need to consider any infections happening after this one
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

        double product = Math.exp(logProduct);

        logP -= Math.log1p(-product);

        return logP;
    }



    private class totalLogProbability implements MultivariateFunction {

        final boolean hasGeography;

        private totalLogProbability(boolean hasGeography){
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

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            CaseToCaseTreeLikelihood c2cTL = (CaseToCaseTreeLikelihood)
                    xo.getChild(CaseToCaseTreeLikelihood.class);
            SpatialKernel kernel = (SpatialKernel) xo.getChild(SpatialKernel.class);
            Parameter transmissionRate = (Parameter) xo.getElementFirstChild(TRANSMISSION_RATE);
            int steps = 50;

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



}
