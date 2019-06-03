package dr.evomodel.coalescent.hmc;

import dr.evomodel.coalescent.GMRFMultilocusSkyrideLikelihood;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;

/**
 * @author Marc A. Suchard
 * @author Mandev Gill
 */
public abstract class GMRFGradient implements GradientWrtParameterProvider {

    private final GMRFMultilocusSkyrideLikelihood skygridLikelihood;
    private final Parameter parameter;

    public GMRFGradient(GMRFMultilocusSkyrideLikelihood skygridLikelihood, Parameter parameter) {
        this.skygridLikelihood = skygridLikelihood;
        this.parameter = parameter;
    }

    @Override
    public Likelihood getLikelihood() {
        return skygridLikelihood;
    }

    @Override
    public Parameter getParameter() {
        return parameter;
    }

    @Override
    public int getDimension() {
        return parameter.getDimension();
    }

    @Override
    abstract public double[] getGradientLogDensity();

    public class WrtLogPopulationSizes extends GMRFGradient {

        public WrtLogPopulationSizes(GMRFMultilocusSkyrideLikelihood skygridLikelihood) {
            super(skygridLikelihood, skygridLikelihood.getPopSizeParameter());
        }

        @Override
        public double[] getGradientLogDensity() {
            return skygridLikelihood.getGradientWrtLogPopulationSize();
        }
    }

    public class WrtPrecision extends GMRFGradient {

        public WrtPrecision(GMRFMultilocusSkyrideLikelihood skygridLikelihood) {
            super(skygridLikelihood, skygridLikelihood.getPrecisionParameter());
        }

        @Override
        public double[] getGradientLogDensity() {
            return skygridLikelihood.getGradientWrtPrecision();
        }
    }

    public class WrtRegressionCoefficients extends GMRFGradient {

        public WrtRegressionCoefficients(GMRFMultilocusSkyrideLikelihood skygridLikelihood) {
            super(skygridLikelihood, skygridLikelihood.getBetaParameter());
        }

        @Override
        public double[] getGradientLogDensity() {
            return skygridLikelihood.getGradientWrtRegressionCoefficients();
        }
    }
}
