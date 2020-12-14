package dr.inference.operators.shrinkage;

import dr.inference.model.Parameter;
import dr.inference.operators.GibbsOperator;
import dr.inference.operators.SimpleMCMCOperator;

import static dr.inferencexml.operators.shrinkage.BayesianBridgeShrinkageOperatorParser.BAYESIAN_BRIDGE_PARSER;

/**
 * @author Marc A. Suchard
 * @author Alex Fisher
 */

public class TemperOperator extends SimpleMCMCOperator implements GibbsOperator {

    private final Parameter parameter;
    private final Parameter target;
    private final double rate;

    private double[] startValues;

    public TemperOperator(Parameter parameter,
                          Parameter target,
                          double rate,
                          double weight) {

        this.parameter = parameter;
        this.target = target;
        this.rate = rate;

        setWeight(weight);
    }

    @Override
    public String getOperatorName() {
        return BAYESIAN_BRIDGE_PARSER;
    }

    @Override
    public double doOperation() {

        if (startValues == null) {
            startValues = parameter.getParameterValues();
        }

        double currentCount = getCount() + 1;

        double objective = currentCount * Math.exp(-rate);

        for (int i = 0; i < parameter.getDimension(); ++i) {

            double currentValue = parameter.getParameterValue(i);
            double targetValue = target.getParameterValue(i % target.getDimension());

            if (objective > 0.1) {
                double startValue = startValues[i];
                double x = startValue + (startValue - targetValue) * Math.exp(-objective);
            } else if (currentValue != targetValue) {
                parameter.setParameterValue(i, targetValue);
            }
        }

        return 0;
    }

}
