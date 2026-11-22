package dr.inference.glm;

import dr.inference.distribution.LogisticRegression;
import dr.inference.model.Parameter;
import dr.math.distributions.SplineBasisMatrix;
import dr.math.distributions.SplineGenerator;

/**
 * Analytic gradient provider for logistic regression log-likelihood with respect to a fixed-effect parameter.
 */
public class LogisticRegressionDesignMatrixGradient extends AbstractLogisticRegressionGradient {

    private final SplineBasisMatrix spline;
    private final SplineGenerator.Derivative generator;
    private final Parameter coefficients;
    private final Parameter delta;

    private double[][] chainRule;

    public LogisticRegressionDesignMatrixGradient(LogisticRegression likelihood, SplineBasisMatrix spline) {
        super(likelihood, spline.getDesignParameter());

        int effectIndex = -1;
        for (int i = 0; i < likelihood.getNumberOfFixedEffects(); ++i) {
            if (likelihood.getDesignMatrix(i) == spline) {
                effectIndex = i;
                break;
            }
        }

        if (effectIndex < 0) {
            throw new IllegalArgumentException("Parameter is not a fixed effect in the supplied logistic regression model");
        }

        this.spline = spline;
        this.generator = new SplineGenerator.Derivative(spline.getDegree(), spline.getIncludeIntercept(),
                spline.getZeroOutOfBound());
        this.coefficients = likelihood.getFixedEffect(effectIndex);
        this.delta = likelihood.getFixedEffectIndicator(effectIndex);
    }

    @Override
    public double[] getGradientLogDensity() {
        final int N = dependent.getDimension();
        final int P = coefficients.getDimension();

        if (chainRule == null) {
            chainRule = new double[N][P];
        }

        final double[] gradient = new double[N];

        double[] xBeta = likelihood.getXBeta();
        generator.fillBasis(chainRule,
                spline.getExpandedKnots(), spline.getDesignParameter(),
                spline.getLowerBoundary(), spline.getUpperBoundary());

        for (int i = 0; i < N; i++) {
            double y = dependent.getParameterValue(i);
            double eta = xBeta[i];
            double derivative = 0.0;
            if (!Double.isNaN(y) && !Double.isNaN(eta)) {
                double p_i = logistic(eta);
                double residual = y - p_i;

                for (int k = 0; k < P; k++) {
                    double increment = chainRule[i][k] * coefficients.getParameterValue(k);
                    if (delta != null) {
                        increment *= delta.getParameterValue(k);
                    }
                    derivative += increment;
                }
                derivative *= residual;
            }
            gradient[i] = derivative;
        }

        return gradient;
    }
}
