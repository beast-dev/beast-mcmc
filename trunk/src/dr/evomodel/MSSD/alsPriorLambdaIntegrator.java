package dr.evomodel.MSSD;

import dr.evolution.alignment.PatternList;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
import dr.inference.loggers.LogColumn;
import dr.inference.model.AbstractModel;
import dr.inference.model.Likelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.math.GammaFunction;
import dr.xml.*;

import java.util.logging.Logger;


/**
 * Package: alsPriorLambdaIntegrator
 * Description:
 * <p/>
 * <p/>
 * Created by
 * Alexander V. Alekseyenko (alexander.alekseyenko@gmail.com)
 * Date: Mar 20, 2008
 * Time: 2:06:55 PM
 */
public class alsPriorLambdaIntegrator extends AbstractModel implements Likelihood {
    protected PatternList patterns;
    protected Parameter deathRate;
    protected Parameter creationRate;
    protected AbstractObservationProcess observationProcess;
    protected BranchRateModel branchRateModel;
    protected int patternCount;
    protected double N, storedN;
    protected double patternWeights[];
    protected boolean patternWeightKnown;
    protected double gammaNorm;

    public alsPriorLambdaIntegrator(BranchRateModel branchRateModel,
                                    AbstractObservationProcess observationProcess,
                                    PatternList patterns) {
        super(MODEL_NAME);
        this.deathRate = observationProcess.getMuParameter();
        this.creationRate = observationProcess.getLamParameter();
        this.observationProcess = observationProcess;
        this.branchRateModel = branchRateModel;
        this.patterns = patterns;
        patternWeightKnown = false;
    }

    public static final String MODEL_NAME = "alsLambdaIntegrator";

    /**
     * The XML parser
     */
    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return MODEL_NAME;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            Logger.getLogger("dr.evolution").info("\n ---------------------------------\nCreating alsLambdaIntegrator model.");
            Logger.getLogger("dr.evolution").info("\tIf you publish results using this prior, please reference:");
            Logger.getLogger("dr.evolution").info("\t\t 1. Ferreira and Suchard (in press) for the conditional reference prior on CTMC scale parameter prior;");
            Logger.getLogger("dr.evolution").info("\t\t 2. Alekseyenko, Lee and Suchard (in submision).");

            PatternList patterns = (PatternList) xo.getChild(PatternList.class);
            BranchRateModel branchRateModel = (BranchRateModel) xo.getChild(BranchRateModel.class);
            if (branchRateModel == null) {
                branchRateModel = new DefaultBranchRateModel();
                Logger.getLogger("dr.evolution").info("\tUsing DefaultBranchRateModel\n---------------------------------\n");
            }
            AbstractObservationProcess observationProcess = (AbstractObservationProcess) xo.getChild(AbstractObservationProcess.class);

            return new alsPriorLambdaIntegrator(branchRateModel, observationProcess, patterns);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element in conjunction with AlternativeSplicingModel integrates the immigration rate parameter.";
        }

        public Class getReturnType() {
            return alsPriorLambdaIntegrator.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(PatternList.class),
                new ElementRule(AbstractObservationProcess.class),
                new ElementRule(BranchRateModel.class, true)
        };
    };

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == patterns) {
            patternWeightKnown = false;
        }
    }

    /**
     * This method is called whenever a parameter is changed.
     * <p/>
     * It is strongly recommended that the model component sets a "dirty" flag and does no
     * further calculations. Recalculation is typically done when the model component is asked for
     * some information that requires them. This mechanism is 'lazy' so that this method
     * can be safely called multiple times with minimal computational cost.
     */
    protected void handleParameterChangedEvent(Parameter parameter, int index) {
        if (parameter == patterns) {
            patternWeightKnown = false;
        }
    }

    /**
     * Additional state information, outside of the sub-model is stored by this call.
     */
    protected void storeState() {
        storedN = N;
    }

    /**
     * After this call the model is guaranteed to have returned its extra state information to
     * the values coinciding with the last storeState call.
     * Sub-models are handled automatically and do not need to be considered in this method.
     */
    protected void restoreState() {
        N = storedN;
    }

    /**
     * This call specifies that the current state is accept. Most models will not need to do anything.
     * Sub-models are handled automatically and do not need to be considered in this method.
     */
    protected void acceptState() {
    }

    /**
     * Get the model.
     *
     * @return the model.
     */
    public Model getModel() {
        return this;  //AUTOGENERATED METHOD IMPLEMENTATION
    }

    /**
     * Get the log likelihood.
     *
     * @return the log likelihood.
     */
    public double getLogLikelihood() {
        return calculateLogLikelihood();  //AUTOGENERATED METHOD IMPLEMENTATION
    }

    public double calculateLogLikelihood() {
        double lam = creationRate.getParameterValue(0);
        double mu = deathRate.getParameterValue(0);
        if (!patternWeightKnown) {
            N = 0.0;
            patternCount = patterns.getPatternCount();
            double patternWeights[] = patterns.getPatternWeights();
            for (int i = 0; i < patternCount; ++i) {
                N += patternWeights[i];
            }
            patternWeightKnown = true;
        }
        double treeWeight = observationProcess.getLogTreeWeight(branchRateModel);
        return GammaFunction.lnGamma(N) - (N - 1) * Math.log(lam) + N * Math.log(mu) - treeWeight - N * Math.log(-treeWeight * mu / lam);
    }

    /**
     * Forces a complete recalculation of the likelihood next time getLikelihood is called
     */
    public void makeDirty() {
        patternWeightKnown = false;
    }

    public LogColumn[] getColumns() {
        return new LogColumn[0];  //AUTOGENERATED METHOD IMPLEMENTATION
    }
}
