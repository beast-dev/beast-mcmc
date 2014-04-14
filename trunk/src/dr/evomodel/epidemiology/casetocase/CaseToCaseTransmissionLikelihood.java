package dr.evomodel.epidemiology.casetocase;

import dr.evomodel.coalescent.DemographicModel;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.model.*;
import dr.math.*;
import dr.xml.*;

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

    private static final boolean DEBUG = false;

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
    private boolean geographyProbKnown;
    private boolean storedGeographyProbKnown;
    private boolean treeProbKnown;
    private boolean storedTreeProbKnown;
    private double logLikelihood;
    private double storedLogLikelihood;
    private double transLogProb;
    private double storedTransLogProb;
    private double normalisation;
    private double storedNormalisation;
    private double geographyLogProb;
    private double storedGeographyLogProb;
    private double treeLogProb;
    private double storedTreeLogProb;
    private final boolean hasGeography;
    private final TransmissionNumericalIntegrator integrator;
    private ArrayList<AbstractCase> sortedCases;
    private ArrayList<AbstractCase> storedSortedCases;
    private TotalProbability totalProbability;

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
        if(transmissionRate.getBounds().getUpperLimit(0)==Double.POSITIVE_INFINITY){
            throw new RuntimeException("Infinite upper limits for transmission rate not implemented yet");
        }
        this.addModel(treeLikelihood);
        this.addVariable(transmissionRate);
        likelihoodKnown = false;
        integrator = new TransmissionNumericalIntegrator(steps);
        hasGeography = spatialKernal!=null;
        totalProbability = new TotalProbability(hasGeography);
        sortCases();
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if(model instanceof CaseToCaseTreeLikelihood){

            treeProbKnown = false;
            if(!(object instanceof AbstractOutbreak) && !(object instanceof DemographicModel)){
                transProbKnown = false;
                normalisationKnown = false;
                geographyProbKnown = false;
                sortedCases = null;
            }
        } else if(model instanceof SpatialKernel){
            transProbKnown = false;

            normalisationKnown = false;

            geographyProbKnown = false;
        }
        likelihoodKnown = false;
    }

    // no need to change the RNG queue unless the normalisation will have changed

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if(variable==transmissionRate){
            transProbKnown = false;

            normalisationKnown = false;
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
        storedGeographyLogProb = geographyLogProb;
        storedGeographyProbKnown = geographyProbKnown;
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
        geographyLogProb = storedGeographyLogProb;
        geographyProbKnown = storedGeographyProbKnown;
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

    public double getLogLikelihood() {

        if(!likelihoodKnown){
            if(!treeProbKnown){
                treeLikelihood.prepareTimings();
            }
            if(!transProbKnown){
                try{
                    double rate = transmissionRate.getParameterValue(0);

                    transLogProb = totalProbability.logEvaluate(rate);
                    integrator.setCentre(rate);
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
                normalisation = integrator.logIntegrate(totalProbability, transmissionRate.getBounds().getLowerLimit(0),
                        transmissionRate.getBounds().getUpperLimit(0));

                normalisationKnown = true;
            }
            if(!geographyProbKnown){
                geographyLogProb = 0;

                for(AbstractCase aCase : outbreak.getCases()){
                    int number = outbreak.getCaseIndex(aCase);

                    double infectionTime = treeLikelihood.getInfectionTime(aCase);
                    AbstractCase parent = treeLikelihood.getInfector(aCase);
                    if(parent!=null){
                        if(treeLikelihood.getInfectiousTime(parent)>infectionTime
                                        || parent.culledYet(infectionTime)) {
                            geographyLogProb += Double.NEGATIVE_INFINITY;
                        } else {
                            double numerator = outbreak.getKernelValue(aCase, parent, spatialKernel);
                            double denominator = 0;

                            for(int i=0; i< outbreak.size(); i++){
                                AbstractCase parentCandidate = outbreak.getCase(i);

                                if(i!=number && treeLikelihood.getInfectiousTime(parentCandidate)<infectionTime
                                        && !parentCandidate.culledYet(infectionTime)){
                                    denominator += (outbreak.getKernelValue(aCase, parentCandidate,
                                            spatialKernel));
                                }
                            }
                            geographyLogProb += Math.log(numerator/denominator);
                        }
                    } else {
                        // probability of first infection given all the timings is 1

                        geographyLogProb += 0;
                    }
                }
                geographyProbKnown = true;
            }

            if(!treeProbKnown){
                treeLogProb = treeLikelihood.getLogLikelihood();
                treeProbKnown = true;
            }
            logLikelihood =  treeLogProb + geographyLogProb + transLogProb - normalisation;
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
        normalisationKnown = false;
        geographyProbKnown = false;
        treeProbKnown = false;
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

    /* todo at the moment we're assuming a uniform prior on the finite bounds of transmissionRate. This should be
    regarded as temporary
     */

    private double caseLogProbability(AbstractCase infectee, double argument){
        double logP = 0;

        if(sortedCases == null){
            sortCases();
        }
        if(treeLikelihood.getInfector(infectee)==null){
            return 0;
        }
        AbstractCase infector = treeLikelihood.getInfector(infectee);

        // probability that the infector infected this case at this time

        double r = argument;

        double[] kernelValues = outbreak.getKernelValues(infectee, spatialKernel);

        if(hasGeography){
            double kernelValue = kernelValues[outbreak.getCaseIndex(infector)];
            r *= kernelValue;
        }

        double infecteeInfected = treeLikelihood.getInfectionTime(infectee);
        double infectorInfectious = treeLikelihood.getInfectiousTime(infector);
        double infecteeExamined = infectee.getExamTime();
        double infectorNoninfectious = infector.getCullTime();

        // need to differentiate between an actual zero probability and a rounding error. This is a zero probability.

        if(infecteeInfected < infectorInfectious || infecteeInfected > infectorNoninfectious){
            throw new BadPartitionException("Illegal partition given known timings");
        }

        logP += Math.log(r);

        // probability that nothing infected this case sooner
        // and

        double sum1 = 0;
        double sum2 = 0;

        for(AbstractCase possibleInfector : sortedCases){


            if(possibleInfector!=infectee){

                double nonInfectorInfectious = treeLikelihood.getInfectiousTime(possibleInfector);

                if(nonInfectorInfectious > infecteeExamined){
                    break;
                }

                double nonInfectorNoninfectious = possibleInfector.getCullTime();

                double rate = hasGeography ? kernelValues[outbreak.getCaseIndex(possibleInfector)] : 1;

                if(nonInfectorInfectious <= infecteeInfected){
                    double lastPossibleInfectionTime = Math.min(nonInfectorNoninfectious, infecteeInfected);
                    sum1 += -rate*(lastPossibleInfectionTime - nonInfectorInfectious);
                }


                double lastPossibleInfectionTime = Math.min(nonInfectorNoninfectious, infecteeExamined);
                sum2 += -rate*(lastPossibleInfectionTime - nonInfectorInfectious);
            }

        }

        // probability that nothing infected this case sooner

        logP += argument*sum1;

        // probability that _something_ infected this case before it was examined

        double product = Math.exp(argument*sum2);

        logP -= Math.log1p(-product);


        // prior probability

        logP -= Math.log(transmissionRate.getBounds().getUpperLimit(0) - transmissionRate.getBounds().getLowerLimit(0));

        return logP;
    }

    private class TotalProbability implements IntegrableUnivariateFunction {

        final boolean hasGeography;

        private TotalProbability(boolean hasGeography){
            this.hasGeography = hasGeography;
        }

        // index 0 is lambda, index 1 if present is alpha

        public double evaluate(double argument) {
            return Math.exp(logEvaluate(argument));
        }

        public double logEvaluate(double argument) {
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

        public double getLowerBound() {
            return 0;
        }

        public double getUpperBound() {
            return transmissionRate.getBounds().getUpperLimit(0);
        }

        public double evaluateIntegral(double a, double b) {
            return integrator.integrate(this, a, b);
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
        ArrayList<LogColumn> columns = new ArrayList<LogColumn>();

        columns.add(new LogColumn.Abstract("trans_LL"){
            protected String getFormattedValue() {
                return String.valueOf(transLogProb - normalisation);
            }
        });

        columns.add(new LogColumn.Abstract("geog_LL"){
            protected String getFormattedValue() {
                return String.valueOf(geographyLogProb);
            }
        });

        columns.addAll(Arrays.asList(treeLikelihood.passColumns()));


        return columns.toArray(new LogColumn[columns.size()]);
    }

    // to prevent underflow:

    public class TransmissionNumericalIntegrator implements Integral{

        private int sampleSize;
        private double centre;

        public TransmissionNumericalIntegrator(int sampleSize){
            this.sampleSize = sampleSize;
        }

        public void setCentre(double value){
            centre = value;
        }

        public double logIntegrate(TotalProbability f, double min, double max){
            if(centre<min || centre>max){
                throw new IllegalArgumentException();
            }

            double logIntegral = Double.NEGATIVE_INFINITY;

            int samplesPerHalf = sampleSize % 2 == 0 ? sampleSize/2 : (sampleSize-1)/2;
            int startOfUpperHalf = sampleSize % 2 == 0 ? sampleSize/2 : (sampleSize+1)/2;

            double bottomHalf = centre - min;
            double bottomStep = bottomHalf/samplesPerHalf;
            double topHalf = max - centre;
            double topStep = topHalf/samplesPerHalf;

            double[] points = new double[sampleSize];

            if(sampleSize % 2 == 1){
                points[(sampleSize - 1)/2]=centre;
            }

            for(int i=0; i<samplesPerHalf; i++){
                points[i]= min + bottomStep/2 + i*bottomStep;
                points[startOfUpperHalf + i]= centre + topStep/2 + i*topStep;
            }

            for (int i = 0; i < sampleSize; i++) {
                logIntegral = LogTricks.logSum(logIntegral, f.logEvaluate(points[i]));
            }


            logIntegral += Math.log((max - min) / (double)sampleSize);
            return logIntegral;
        }

        public double integrate(UnivariateFunction f, double min, double max) {
            return Math.exp(logIntegrate((TotalProbability)f, min, max));
        }
    }

}
