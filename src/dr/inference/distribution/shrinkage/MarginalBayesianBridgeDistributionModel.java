package dr.inference.distribution.shrinkage;

import dr.inference.model.Parameter;
import dr.math.distributions.MarginalizedAlphaStableDistribution;

/**
 * @author Marc A. Suchard
 * @author Akihiko Nishimura
 */

public class MarginalBayesianBridgeDistributionModel extends BayesianBridgeDistributionModel {

    public MarginalBayesianBridgeDistributionModel(Parameter globalScale,
                                                   Parameter exponent,
                                                   int dim) {
        super(globalScale, exponent, dim);
    }

    @Override
    public Parameter getLocalScale() { return null; }

    @Override
    public Parameter getSlabWidth() { return null; }

    @Override
    double[] gradientLogPdf(double[] x) {
        final int dim = x.length;
        final double scale = globalScale.getParameterValue(0);
        final double alpha = exponent.getParameterValue(0);

        double[] gradient = new double[dim];
        for (int i = 0; i < dim; ++i) {
            gradient[i] = MarginalizedAlphaStableDistribution.gradLogPdf(x[i], scale, alpha);
        }
        return gradient;
    }

    @Override
    public double logPdf(double[] v) {
        final double scale = globalScale.getParameterValue(0);
        final double alpha = exponent.getParameterValue(0);

        double sum = 0.0;
        for (double x : v) {
            sum += MarginalizedAlphaStableDistribution.logPdf(x, scale, alpha);
        }
        return sum;
    }
}