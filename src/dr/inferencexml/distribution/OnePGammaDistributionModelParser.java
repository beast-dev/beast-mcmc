package dr.inferencexml.distribution;

import dr.inference.distribution.GammaDistributionModel;
import dr.inference.distribution.ParametricDistributionModel;
import dr.inference.model.Parameter;

/**
 */
public class OnePGammaDistributionModelParser extends DistributionModelParser {
    
    public String getParserName() {
        return GammaDistributionModel.ONE_P_GAMMA_DISTRIBUTION_MODEL;
    }

    ParametricDistributionModel parseDistributionModel(Parameter[] parameters, double offset) {
        return new GammaDistributionModel(parameters[0]);
    }

    public String[] getParameterNames() {
        return new String[]{SHAPE};
    }

    public String getParserDescription() {
        return "A model of a one parameter gamma distribution.";
    }

    public boolean allowOffset() {
        return false;
    }

    public Class getReturnType() {
        return GammaDistributionModel.class;
    }
}
