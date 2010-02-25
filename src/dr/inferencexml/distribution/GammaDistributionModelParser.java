package dr.inferencexml.distribution;

import dr.inference.distribution.GammaDistributionModel;
import dr.inference.distribution.ParametricDistributionModel;
import dr.inference.model.Parameter;

/**
 */
public class GammaDistributionModelParser extends DistributionModelParser {

    public String getParserName() {
        return GammaDistributionModel.GAMMA_DISTRIBUTION_MODEL;
    }

    ParametricDistributionModel parseDistributionModel(Parameter[] parameters, double offset) {
        return new GammaDistributionModel(parameters[0], parameters[1]);
    }

    public String[] getParameterNames() {
        return new String[]{SHAPE, SCALE};
    }

    public String getParserDescription() {
        return "A model of a gamma distribution.";
    }

    public boolean allowOffset() {
        return false;
    }

    public Class getReturnType() {
        return GammaDistributionModel.class;
    }
}
