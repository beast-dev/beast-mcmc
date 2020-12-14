package dr.math.distributions;

import dr.math.UnivariateFunction;

public class ParetoDistribution implements Distribution {

    private double scale;
    private double shape;

    public ParetoDistribution(double scale, double shape) {
        if (scale <= 0 || shape <= 0) {
            throw new RuntimeException("Shape and scale must be positive.");
        }
        this.scale = scale;
        this.shape = shape;
    }

    @Override
    public double pdf(double x) {
        return x > scale ? shape * Math.pow(scale, shape) / Math.pow(x, shape + 1) : 0.0;
    }

    @Override
    public double logPdf(double x) {
        if (x < scale) return Double.NEGATIVE_INFINITY;
        return Math.log(shape) + shape * Math.log(scale) - (shape + 1) * Math.log(x);
    }

    @Override
    public double cdf(double x) {
        if (x < scale) return 0.0;
        return 1 - Math.pow(scale / x, shape);
    }

    @Override
    public double quantile(double p) {
        return 1.0 / Math.pow(1.0 - p, 1 / shape);
    }//todo

    @Override
    public double mean() {
        if (shape <= 2) return Double.POSITIVE_INFINITY;
        else return scale * shape / (shape - 1);
    }

    @Override
    public double variance() {
        if (shape <= 2) return Double.POSITIVE_INFINITY;
        else return scale * scale * shape / (shape - 1) * (shape - 1) * (shape - 2);
    }

    @Override
    public UnivariateFunction getProbabilityDensityFunction() {
        return null;
    }
}
