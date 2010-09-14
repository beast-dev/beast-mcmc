package dr.inferencexml.distribution;

import dr.inference.distribution.ExponentialDistributionModel;
import dr.inference.distribution.ParametricDistributionModel;
import dr.inference.model.Parameter;

/**
 */
public class ExponentialDistributionModelParser extends DistributionModelParser {

    public String getParserName() {
        return ExponentialDistributionModel.EXPONENTIAL_DISTRIBUTION_MODEL;
    }

    ParametricDistributionModel parseDistributionModel(Parameter[] parameters, double offset) {
        return new ExponentialDistributionModel(parameters[0], offset);
    }

    public String[] getParameterNames() {
        return new String[]{MEAN};
    }

    public String getParserDescription() {
        return "A model of an exponential distribution.";
    }

    public boolean allowOffset() {
        return true;
    }

    public Class getReturnType() {
        return ExponentialDistributionModel.class;
    }
}
