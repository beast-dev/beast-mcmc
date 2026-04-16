package dr.inference.glm;

import dr.inference.distribution.LogisticRegression;
import dr.inference.model.DesignMatrix;
import dr.inference.model.Parameter;

/**
 * Analytic gradient provider for logistic regression log-likelihood with respect to a fixed-effect parameter.
 */
public class LogisticRegressionGradientWrtParameter extends AbstractLogisticRegressionGradient {

    private final DesignMatrix designMatrix;
    private final Parameter delta;

    public LogisticRegressionGradientWrtParameter(LogisticRegression likelihood, Parameter parameter) {
        super(likelihood, parameter);

        int effectIndex = likelihood.getEffectNumber(parameter);
        if (effectIndex < 0) {
            throw new IllegalArgumentException("Parameter is not a fixed effect in the supplied logistic regression model");
        }

        this.designMatrix = likelihood.getDesignMatrix(effectIndex);
        this.delta = likelihood.getFixedEffectIndicator(effectIndex);

    }

    @Override
    public double[] getGradientLogDensity() {
        final int n = dependent.getDimension();
        final int p = parameter.getDimension();
        final double[] gradient = new double[p];

        double[] xBeta = likelihood.getXBeta();

        for (int i = 0; i < n; i++) {
            double y = dependent.getParameterValue(i);
            if (!Double.isNaN(y)) {
                double eta = xBeta[i];
                double p_i = logistic(eta);
                double residual = y - p_i;

                for (int k = 0; k < p; k++) {
                    double multiplier = designMatrix.getParameterValue(i, k);
                    if (delta != null) {
                        multiplier *= delta.getParameterValue(k);
                    }
                    gradient[k] += multiplier * residual;
                }
            }
        }

        return gradient;
    }
}
