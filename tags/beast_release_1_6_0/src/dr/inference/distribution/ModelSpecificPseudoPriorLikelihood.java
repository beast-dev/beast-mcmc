package dr.inference.distribution;

import dr.math.distributions.Distribution;
import dr.inference.model.Parameter;

/**
 * @author Chieh-Hsi Wu
 */
public class ModelSpecificPseudoPriorLikelihood extends DistributionLikelihood{

    private int[] models;
    protected Distribution prior;
    protected Distribution pseudoPrior;
    private Parameter modelIndicator;

    public ModelSpecificPseudoPriorLikelihood(
            Distribution prior,
            Distribution pseudoPrior,
            Parameter modelIndicator,
            int[] models){
        super(prior);
        this.prior = prior;
        this.pseudoPrior = pseudoPrior;
        this.models = models;
        this.modelIndicator = modelIndicator;
    }


   /**
     * Calculate the log likelihood of the current state.
     *
     * @return the log likelihood.
     */
    public double calculateLogLikelihood() {
        boolean inModel = false;
        int modelCode = (int)modelIndicator.getParameterValue(0);
        for(int i = 0; i < models.length; i++){
            if(models[i] == modelCode){
                inModel = true;
                break;
            }
        }

        if(inModel){
            distribution = prior;
        }else{
            distribution = pseudoPrior;
        }
        double logL = super.calculateLogLikelihood();
        return logL;
    }

}
