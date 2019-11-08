package dr.inference.distribution.shrinkage;

import dr.inference.model.Parameter;
import dr.math.distributions.MarginalizedAlphaStableDistribution;

public class OldMarginalBayesianBridge extends OldBayesianBridgeLikelihood {

    public OldMarginalBayesianBridge(Parameter coefficients,
                                     Parameter globalScale,
                                     Parameter exponent) {
        super(coefficients, globalScale, exponent);
    }

    double calculateLogLikelihood() {

        final double scale = globalScale.getParameterValue(0);
        final double alpha = exponent.getParameterValue(0);

        double sum = 0.0;
        for (int i = 0; i < dim; ++i) {
            sum += MarginalizedAlphaStableDistribution.logPdf(coefficients.getParameterValue(i), scale, alpha);
        }
        return sum;
    }

    @Override
    double[] calculateGradientLogDensity() {

        final double scale = globalScale.getParameterValue(0);
        final double alpha = exponent.getParameterValue(0);

        double[] gradient = new double[dim];
        for (int i = 0; i < dim; ++i) {
            gradient[i] = MarginalizedAlphaStableDistribution.gradLogPdf(coefficients.getParameterValue(i),
                    scale, alpha);
        }
        return gradient;
    }

    @Override
    public Parameter getLocalScale() { return null; }
}