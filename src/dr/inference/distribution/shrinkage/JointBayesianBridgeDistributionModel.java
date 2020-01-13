package dr.inference.distribution.shrinkage;

import dr.inference.model.Parameter;
import dr.math.distributions.NormalDistribution;

/**
 * @author Marc A. Suchard
 * @author Akihiko Nishimura
 */

public class JointBayesianBridgeDistributionModel extends BayesianBridgeDistributionModel {

    public JointBayesianBridgeDistributionModel(Parameter globalScale,
                                                Parameter localScale,
                                                Parameter exponent,
                                                Parameter slabWidth,
                                                int dim) {
        super(globalScale, exponent, dim);
        this.localScale = localScale;
        this.slabWidth = slabWidth;

        if (dim != localScale.getDimension()) {
            throw new IllegalArgumentException("Invalid dimensions");
        }

        addVariable(localScale);
    }

    @Override
    public Parameter getLocalScale() { return localScale; }

    @Override
    public Parameter getSlabWidth() { return slabWidth; }

    @Override
    double[] gradientLogPdf(double[] x) {

        double[] gradient = new double[dim];
        for (int i = 0; i < dim; ++i) {
            gradient[i] = NormalDistribution.gradLogPdf(x[i], 0, getStandardDeviation(i));
        }
        return gradient;
    }

    @Override
    public double logPdf(double[] x) {

        double pdf = 0.0;
        for (int i = 0; i < dim; ++i) {
            pdf += NormalDistribution.logPdf(x[i], 0, getStandardDeviation(i));
        }

        // TODO Add density of localScale variables

        return pdf;
    }

    private double getStandardDeviation(int index) {
        double globalLocalProduct = globalScale.getParameterValue(0) * localScale.getParameterValue(index);
        if (slabWidth != null) {
            double ratio = globalLocalProduct / slabWidth.getParameterValue(0);
            globalLocalProduct /= Math.sqrt(1.0 + ratio * ratio);
        }
        return globalLocalProduct;
    }

    private final Parameter localScale;
    private final Parameter slabWidth;
}