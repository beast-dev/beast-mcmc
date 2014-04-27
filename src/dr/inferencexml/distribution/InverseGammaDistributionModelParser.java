package dr.inferencexml.distribution;

import dr.inference.distribution.GammaDistributionModel;
import dr.inference.distribution.InverseGammaDistributionModel;
import dr.inference.distribution.ParametricDistributionModel;
import dr.inference.model.Parameter;

/**
 */
public class InverseGammaDistributionModelParser extends DistributionModelParser {

    public String getParserName() {
        return InverseGammaDistributionModel.INVERSE_GAMMA_DISTRIBUTION_MODEL;
    }

    ParametricDistributionModel parseDistributionModel(Parameter[] parameters, double offset) {
        return new InverseGammaDistributionModel(parameters[0], parameters[1]);
    }

    public String[] getParameterNames() {
        return new String[]{SHAPE, SCALE};
    }

    public String getParserDescription() {
        return "A model of an inverese gamma distribution.";
    }

    public boolean allowOffset() {
        return false;
    }

    public Class getReturnType() {
        return InverseGammaDistributionModel.class;
    }
}
