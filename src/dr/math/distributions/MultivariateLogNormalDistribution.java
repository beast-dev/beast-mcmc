package dr.math.distributions;

import dr.inference.model.GradientProvider;

public class MultivariateLogNormalDistribution implements MultivariateDistribution, GradientProvider {

    public static final String TYPE = "MultivariateLogNormal";
    private final LogNormalDistribution logNormalDistribution;

    private final int dim;

    public MultivariateLogNormalDistribution(LogNormalDistribution logNormalDistribution, int dim) {
        this.logNormalDistribution = logNormalDistribution;
        this.dim = dim;
    }

    @Override
    public int getDimension() {
        return dim;
    }

    @Override
    public double[] getGradientLogDensity(Object x) {
        return logNormalDistribution.getGradientLogDensity(x);
    }

    @Override
    public double logPdf(double[] x) {
        double logPdf = 0;
        for (int i = 0; i < dim; i++) {
            logPdf += logNormalDistribution.logPdf(x[i]);
        }
        return logPdf;
    }

    @Override
    public double[][] getScaleMatrix() {
        return new double[0][];
    }

    @Override
    public double[] getMean() {
        double[] meanVector = new double[dim];
        for (int i = 0; i < dim; i++) {
            meanVector[i] = logNormalDistribution.M;
        }
        return meanVector;
    }

    @Override
    public String getType() {
        return null;
    }
}
