package dr.inference.prior;

import dr.evolution.alignment.PatternList;
import dr.inference.loggers.LogColumn;
import dr.inference.model.AbstractModel;
import dr.inference.model.Likelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.math.GammaFunction;
import dr.xml.*;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.branchratemodel.BranchRateModel;

import java.util.logging.Logger;


/**
 * Package: ImmigrationIntegrationForAlternativeSplicingPrior
 * Description:
 * <p/>
 * <p/>
 * Created by
 * Alexander V. Alekseyenko (alexander.alekseyenko@gmail.com)
 * Date: Mar 20, 2008
 * Time: 2:06:55 PM
 */
public class ImmigrationIntegrationForAlternativeSplicingPrior extends AbstractModel implements Likelihood {
    protected PatternList patterns;
    protected Parameter deathRate;
    protected Parameter creationRate;
    protected int patternCount;
    protected double N, storedN;
    protected double patternWeights[];
    protected boolean patternWeightKnown;
    protected double gammaNorm;

    public ImmigrationIntegrationForAlternativeSplicingPrior(Parameter deathRate,
                                                             Parameter creationRate, PatternList patterns) {
        super(MODEL_NAME);
        this.deathRate = deathRate;
        this.creationRate = creationRate;
        this.patterns = patterns;
        patternWeightKnown=false;
    }


    public static final String MODEL_NAME="immigrationIntegrationPrior";
    public static final String DEATHPARAMETER="deathRate";
    public static final String CREATIONPARAMETER = "immigrationRate";

    /**
	 * The XML parser
	 */
	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

		public String getParserName() { return MODEL_NAME; }

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            PatternList patterns= (PatternList) xo.getChild(PatternList.class);
            Parameter creationRate = (Parameter)xo.getSocketChild(CREATIONPARAMETER);
            Parameter deathRate = (Parameter)xo.getSocketChild(DEATHPARAMETER);

            Logger.getLogger("dr.evolution").info("\n ---------------------------------\nCreating ImmigrationIntegrationPrior model.");
            Logger.getLogger("dr.evolution").info("\tIf you publish results using this prior, please reference:");
            Logger.getLogger("dr.evolution").info("\t\t 1. Ferreira and Suchard (in press) for the conditional reference prior on CTMC scale parameter prior;");
            Logger.getLogger("dr.evolution").info("\t\t 2. Alekseyenko, Lee and Suchard (in submision).\n---------------------------------\n");

            return new ImmigrationIntegrationForAlternativeSplicingPrior(deathRate,creationRate,patterns);
        }


        //************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************

		public String getParserDescription() {
			return "This element in conjunction with AlternativeSplicingModel integrates the immigration rate parameter.";
		}

		public Class getReturnType() { return ImmigrationIntegrationForAlternativeSplicingPrior.class; }

		public XMLSyntaxRule[] getSyntaxRules() { return rules; }

		private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
            new ElementRule(PatternList.class),
            new ElementRule(DEATHPARAMETER, new XMLSyntaxRule[] { new ElementRule(Parameter.class) }),
            new ElementRule(CREATIONPARAMETER, new XMLSyntaxRule[] { new ElementRule(Parameter.class) }),
        };
	};

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if(model==patterns){
            patternWeightKnown=false;
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
        if(parameter==patterns){
            patternWeightKnown=false;
        }
    }

    /**
     * Additional state information, outside of the sub-model is stored by this call.
     */
    protected void storeState() {
          storedN=N;
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

    public double calculateLogLikelihood(){
        double lam=creationRate.getParameterValue(0);
        double mu=deathRate.getParameterValue(0);
        if(!patternWeightKnown){
            N=0.0;
            patternCount = patterns.getPatternCount();
            double patternWeights[]=patterns.getPatternWeights();
            for(int i=0; i< patternCount; ++i){
                N+=patternWeights[i];
            }
            patternWeightKnown=true;
        }

        return GammaFunction.lnGamma(N)-(N-1)*Math.log(lam)+N*Math.log(mu)+lam/mu;
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
