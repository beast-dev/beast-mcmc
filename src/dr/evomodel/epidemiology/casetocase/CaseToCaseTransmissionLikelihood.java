package dr.evomodel.epidemiology.casetocase;

import dr.evomodel.coalescent.DemographicModel;
import dr.inference.distribution.GammaDistributionModel;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.model.*;
import dr.math.*;
import dr.math.distributions.GammaDistribution;
import dr.xml.*;

import java.math.BigDecimal;
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
    private Parameter transmissionRate;
    private GammaDistributionModel transmissionRatePrior;
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

    private double betaGammaThing;
    private double storedBetaGammaThing;

    private final boolean hasGeography;
    private ArrayList<AbstractCase> sortedCases;
    private ArrayList<AbstractCase> storedSortedCases;
//    private F f;

    public static final String CASE_TO_CASE_TRANSMISSION_LIKELIHOOD = "caseToCaseTransmissionLikelihood";

    public CaseToCaseTransmissionLikelihood(String name, AbstractOutbreak outbreak,
                                            CaseToCaseTreeLikelihood treeLikelihood, SpatialKernel spatialKernal,
                                            Parameter transmissionRate, GammaDistributionModel transmissionRatePrior){
        super(name);
        this.outbreak = outbreak;
        this.treeLikelihood = treeLikelihood;
        this.spatialKernel = spatialKernal;
        if(spatialKernal!=null){
            this.addModel(spatialKernal);
        }
        this.transmissionRate = transmissionRate;
        this.transmissionRatePrior = transmissionRatePrior;
        this.addModel(treeLikelihood);
        this.addVariable(transmissionRate);
        likelihoodKnown = false;
        hasGeography = spatialKernal!=null;
//        f = new F(hasGeography);
        sortCases();
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if(model instanceof CaseToCaseTreeLikelihood){

            treeProbKnown = false;
            if(!(object instanceof DemographicModel)){
                transProbKnown = false;
                normalisationKnown = false;
                geographyProbKnown = false;
                sortedCases = null;
            }
        } else if(model instanceof SpatialKernel){

            transProbKnown = false;
            normalisationKnown = false;
            geographyProbKnown = false;

        } else if(model instanceof AbstractOutbreak){

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
        storedBetaGammaThing = betaGammaThing;
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
        betaGammaThing = storedBetaGammaThing;
        sortedCases = storedSortedCases;
    }

    protected void acceptState() {
        // nothing to do
    }

    public SpatialKernel getSpatialKernel(){
        return spatialKernel;
    }

    public Model getModel() {
        return this;
    }

    public CaseToCaseTreeLikelihood getTreeLikelihood(){
        return treeLikelihood;
    }

    public double getLogLikelihood() {

        if(DEBUG){
            treeLikelihood.debugOutputTree("blah.nex", true);
        }

        if(!likelihoodKnown) {
            if (!treeProbKnown) {
                treeLikelihood.prepareTimings();
            }
            if (!transProbKnown) {
                double rate = transmissionRate.getParameterValue(0);

                try {
                    double K = getK();

                    // not necessary to actually add it in because it cancels, but need to check the exception
                    getLogD();
                    double E = getE();
                    //double logF = f.logEvaluate(rate);

                    transLogProb = K * Math.log(rate) - E * rate; // + logF;

                    transProbKnown = true;

                    if (!normalisationKnown) {

                        normalisation = GammaFunction.lnGamma(K + 1) - (K+1)*Math.log(E);

                        betaGammaThing = (K+1)/E;

                        //integrator.setAlphaAndB(K,E);

                        //normalisation = integrator.logIntegrate(f, transmissionRate.getBounds().getLowerLimit(0));

                        // not necessary because it cancels
                        //normalisation += logD;

                        normalisationKnown = true;
                    }
                } catch (BadPartitionException e) {
                    transLogProb = Double.NEGATIVE_INFINITY;
                    logLikelihood = Double.NEGATIVE_INFINITY;
                    likelihoodKnown = true;
                    return logLikelihood;
                }


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

            // just reject states where these round to +INF

            if(transLogProb == Double.POSITIVE_INFINITY){
                System.out.println("TransLogProb +INF");
                return Double.NEGATIVE_INFINITY;
            }
            if(geographyLogProb == Double.POSITIVE_INFINITY){
                System.out.println("GeogLogProb +INF");
                return Double.NEGATIVE_INFINITY;
            }
            if(normalisation == Double.NEGATIVE_INFINITY){
                System.out.println("Normalisation +INF");
                return Double.NEGATIVE_INFINITY;
            }
            if(treeLogProb == Double.POSITIVE_INFINITY){
                System.out.println("TreeLogProb +INF");
                return Double.NEGATIVE_INFINITY;
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

    private double getLogD(){

        if(sortedCases == null){
            sortCases();
        }
        double logD = 0;

        for (AbstractCase infectee : sortedCases) {
            double infecteeInfected = treeLikelihood.getInfectionTime(infectee);
            double infecteeInfectious = treeLikelihood.getInfectiousTime(infectee);
            double infecteeNoninfectious = infectee.getCullTime();

            if(infecteeInfected > infecteeInfectious || infecteeInfectious > infecteeNoninfectious){
                throw new BadPartitionException("Illegal partition given known timings");
            }

            AbstractCase infector = treeLikelihood.getInfector(infectee);
            if (infector != null) {
                double infectorInfectious = treeLikelihood.getInfectiousTime(infector);
                double infectorNoninfectious = infector.getCullTime();

                if(infecteeInfected < infectorInfectious || infecteeInfected > infectorNoninfectious){
                    throw new BadPartitionException("Illegal partition given known timings");
                }

                if(hasGeography) {
                    logD += Math.log(outbreak.getKernelValue(infectee, infector, spatialKernel));
                }
            }
        }


        if(transmissionRatePrior!=null) {
            logD += -transmissionRatePrior.getShape() * Math.log(transmissionRatePrior.getScale());
            logD += -GammaFunction.lnGamma(transmissionRatePrior.getShape());
        }

        return logD;

    }

    private double getD(){
        return Math.exp(getLogD());
    }

    private double getE(){
        double E = 0;

        if(sortedCases == null){
            sortCases();
        }
        for(AbstractCase infectee : sortedCases){
            AbstractCase infector = treeLikelihood.getInfector(infectee);
            if (infector != null) {

                double[] kernelValues = outbreak.getKernelValues(infectee, spatialKernel);

                double infecteeInfected = treeLikelihood.getInfectionTime(infectee);

                for (AbstractCase possibleInfector : sortedCases) {

                    if (possibleInfector != infectee) {

                        double nonInfectorInfected = treeLikelihood.getInfectionTime(possibleInfector);
                        double nonInfectorInfectious = treeLikelihood.getInfectiousTime(possibleInfector);
                        double nonInfectorNoninfectious = possibleInfector.getCullTime();

                        if (nonInfectorInfected > infecteeInfected) {
                            break;
                        }

                        double kernelValue = hasGeography ? kernelValues[outbreak.getCaseIndex(possibleInfector)] : 1;

                        if (nonInfectorInfectious <= infecteeInfected) {
                            double lastPossibleInfectionTime = Math.min(nonInfectorNoninfectious, infecteeInfected);
                            E += kernelValue * (lastPossibleInfectionTime - nonInfectorInfectious);
                        }

                    }
                }
            }
        }

        if(transmissionRatePrior!=null) {
            E += 1 / (transmissionRatePrior.getScale());
        }

        return E;

    }

    public double getK(){
        if(transmissionRatePrior != null) {
            return (transmissionRatePrior.getShape() - 1) + outbreak.size()-1;
        } else {
            return outbreak.size()-1;
        }
    }

//    private class F extends UnivariateFunction.AbstractLogEvaluatableUnivariateFunction {
//
//        final boolean hasGeography;
//
//        private F(boolean hasGeography){
//            this.hasGeography = hasGeography;
//        }
//
//        // index 0 is lambda, index 1 if present is alpha
//
//        public double evaluate(double argument) {
//            return Math.exp(logEvaluate(argument));
//        }
//
//        public double logEvaluate(double argument) {
//            if(sortedCases == null){
//                sortCases();
//            }
//            double logF = 0;
//
//            for(AbstractCase infectee : sortedCases){
//
//                AbstractCase infector = treeLikelihood.getInfector(infectee);
//
//                if(infector != null) {
//                    double sum = 0;
//
//                    double[] kernelValues = outbreak.getKernelValues(infectee, spatialKernel);
//
//                    double infecteeExamined = infectee.getExamTime();
//
//                    for (AbstractCase possibleInfector : sortedCases) {
//                        if (possibleInfector != infectee) {
//
//                            double nonInfectorInfected = treeLikelihood.getInfectionTime(possibleInfector);
//                            double nonInfectorInfectious = treeLikelihood.getInfectiousTime(possibleInfector);
//                            double nonInfectorNoninfectious = possibleInfector.getCullTime();
//
//                            if (nonInfectorInfected > infecteeExamined) {
//                                break;
//                            }
//
//                            double rate = hasGeography ? kernelValues[outbreak.getCaseIndex(possibleInfector)] : 1;
//
//                            if (nonInfectorInfectious <= infecteeExamined) {
//                                double lastPossibleInfectionTime = Math.min(nonInfectorNoninfectious, infecteeExamined);
//
//                                sum += -rate * (lastPossibleInfectionTime - nonInfectorInfectious);
//                            }
//                        }
//
//                    }
//                    sum *= argument;
//
//                    double normExp = Math.exp(sum);
//
//                    double logTerm;
//
//                    if(normExp!=1){
//                        logTerm = Math.log1p(-normExp);
//                    } else {
//                        try {
//                            logTerm = handleDenominatorUnderflow(sum);
//                        } catch(IllegalArgumentException e){
//                            throw new RuntimeException("HandleDenominatorUnderflow failed, input = "+sum);
//                        }
//                    }
//
//                    logF += logTerm;
//                }
//            }
//            return (outbreak.size()-1)*Math.log(argument)-logF;
//        }
//
//
//
//        public int getNumArguments() {
//            return 1;
//        }
//
//        public double getLowerBound() {
//            return 0;
//        }
//
//        public double getUpperBound() {
//            return transmissionRate.getBounds().getUpperLimit(0);
//        }
//
//        public double evaluateIntegral(double a, double b) {
//            return integrator.integrate(this, a, b);
//        }
//    }

    private static double handleDenominatorUnderflow(double input){
        BigDecimal bigDec = new BigDecimal(input);
        BigDecimal expBigDec = BigDecimalUtils.exp(bigDec, bigDec.scale());
        BigDecimal one = new BigDecimal(1.0);
        BigDecimal oneMinusExpBigDec = one.subtract(expBigDec);
        BigDecimal logOneMinusExpBigDec = BigDecimalUtils.ln(oneMinusExpBigDec, oneMinusExpBigDec.scale());
        return logOneMinusExpBigDec.doubleValue();
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public static final String TRANSMISSION_RATE = "transmissionRate";
        public static final String INTEGRATOR_STEPS = "integratorSteps";
        public static final String TRANSMISSION_RATE_PRIOR = "transmissionRatePrior";

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            CaseToCaseTreeLikelihood c2cTL = (CaseToCaseTreeLikelihood)
                    xo.getChild(CaseToCaseTreeLikelihood.class);
            SpatialKernel kernel = (SpatialKernel) xo.getChild(SpatialKernel.class);
            Parameter transmissionRate = (Parameter) xo.getElementFirstChild(TRANSMISSION_RATE);

            GammaDistributionModel transmissionRatePrior = null;

            if(xo.hasChildNamed(TRANSMISSION_RATE_PRIOR)) {
                transmissionRatePrior = (GammaDistributionModel) xo.getElementFirstChild(TRANSMISSION_RATE_PRIOR);
            }


            return new CaseToCaseTransmissionLikelihood(CASE_TO_CASE_TRANSMISSION_LIKELIHOOD, c2cTL.getOutbreak(),
                    c2cTL, kernel, transmissionRate, transmissionRatePrior);
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
                new ElementRule(TRANSMISSION_RATE_PRIOR, GammaDistributionModel.class, "A gamma prior on the base" +
                        "transmission rate", true)
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

        columns.add(new LogColumn.Abstract("betaGammaThing"){
            protected String getFormattedValue() {
                return String.valueOf(betaGammaThing);
            }
        });

        columns.addAll(Arrays.asList(treeLikelihood.passColumns()));


        return columns.toArray(new LogColumn[columns.size()]);
    }



}
