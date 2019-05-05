package dr.inference.distribution.shrinkage;

import dr.inference.model.Parameter;
import dr.math.distributions.NormalDistribution;

public class JointBayesianBridgeDistributionModel extends BayesianBridgeDistributionModel {

    public JointBayesianBridgeDistributionModel(Parameter globalScale,
                                                Parameter localScale,
                                                Parameter exponent) {
        super(globalScale, exponent);
        this.localScale = localScale;
        this.dim = localScale.getDimension();

        addVariable(localScale);
    }

    @Override
    public Parameter getLocalScale() { return localScale; }

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
        return globalScale.getParameterValue(0) * localScale.getParameterValue(index);
    }

    private final Parameter localScale;
    private final int dim;
}