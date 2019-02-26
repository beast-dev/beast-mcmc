package dr.inference.distribution.shrinkage;

import dr.inference.model.Parameter;
import dr.math.distributions.NormalDistribution;

public class JointBayesianBridge extends BayesianBridgeLikelihood {

    public JointBayesianBridge(Parameter coefficients,
                               Parameter globalScale,
                               Parameter localScale,
                               Parameter exponent) {
        super(coefficients, globalScale, exponent);
        this.localScale = localScale;

        addVariable(localScale);
    }

    double calculateLogLikelihood() {
        double sum = 0.0;
        for (int i = 0; i < dim; ++i) {
            sum += NormalDistribution.logPdf(coefficients.getParameterValue(i), 0, getStandardDeviation(i));
        }
        return sum;
    }

    @Override
    double[] calculateGradientLogDensity() {

        double[] gradient = new double[dim];
        for (int i = 0; i < dim; ++i) {
            gradient[i] = NormalDistribution.gradLogPdf(coefficients.getParameterValue(i),
                    0, getStandardDeviation(i));
        }
        return gradient;
    }

    @Override
    public Parameter getLocalScale() { return localScale; }

    private double getStandardDeviation(int index) {
        return globalScale.getParameterValue(0) * localScale.getParameterValue(index);
    }

    private final Parameter localScale;
}