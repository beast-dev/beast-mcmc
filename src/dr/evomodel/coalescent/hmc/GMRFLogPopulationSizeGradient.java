package dr.evomodel.coalescent.hmc;

import dr.evomodel.coalescent.GMRFMultilocusSkyrideLikelihood;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;

/**
 * @author Marc A. Suchard
 */
public class GMRFLogPopulationSizeGradient implements GradientWrtParameterProvider {

    private final GMRFMultilocusSkyrideLikelihood skygridLikelihood;
    private final Parameter logPopParameter;

    public GMRFLogPopulationSizeGradient(GMRFMultilocusSkyrideLikelihood skygridLikelihood) {
        this.skygridLikelihood = skygridLikelihood;
        this.logPopParameter = skygridLikelihood.getPopSizeParameter();
    }

    @Override
    public Likelihood getLikelihood() {
        return skygridLikelihood;
    }

    @Override
    public Parameter getParameter() {
        return logPopParameter;
    }

    @Override
    public int getDimension() {
        return logPopParameter.getDimension();
    }

    @Override
    public double[] getGradientLogDensity() {
        return skygridLikelihood.getGradientWrtLogPopulationSize();
    }
}
