package dr.inference.glm;

import dr.inference.distribution.LogisticRegression;
import dr.inference.model.Parameter;

/**
 * Analytic gradient provider for logistic regression log-likelihood with respect to random effects.
 * For each observation i, the random effect appears linearly in the linear predictor:
 * eta_i = X_i * beta + randomEffect_i
 * The gradient for a random effect is simply: dL/d(randomEffect_i) = y_i - p_i
 */
public class LogisticRegressionRandomEffectsGradientWrtParameter extends AbstractLogisticRegressionGradient {

    public LogisticRegressionRandomEffectsGradientWrtParameter(LogisticRegression likelihood, Parameter parameter) {
        super(likelihood, parameter);

        boolean found = false;
        for (int i = 0; i < likelihood.getNumberOfRandomEffects(); i++) {
            if (likelihood.getRandomEffect(i) == parameter) {
                found = true;
                break;
            }
        }
        if (!found) {
            throw new IllegalArgumentException("Parameter is not a random effect in the supplied logistic regression model");
        }

    }

    @Override
    public double[] getGradientLogDensity() {
        final int n = dependent.getDimension();
        final double[] gradient = new double[n];

        double[] xBeta = likelihood.getXBeta();

        for (int i = 0; i < n; i++) {
            double y = dependent.getParameterValue(i);
            if (!Double.isNaN(y)) {
                double eta = xBeta[i];
                double p_i = logistic(eta);
                gradient[i] = y - p_i;
            }
        }

        return gradient;
    }
}
