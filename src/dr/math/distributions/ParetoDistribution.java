package dr.math.distributions;

import dr.math.UnivariateFunction;

public class ParetoDistribution implements Distribution{
    @Override
    public double pdf(double x) {
        return x > 1 ? 1 / x * x : 0.0;
    }

    @Override
    public double logPdf(double x) {
        if (x < 1) return Double.NEGATIVE_INFINITY;
        return -2 * Math.log(x);
    }

    @Override
    public double cdf(double x) { //todo: finish
        return 0;
    }

    @Override
    public double quantile(double y) {
        return 1.0 / (1.0 - y);
    }

    @Override
    public double mean() {
        return Double.POSITIVE_INFINITY;
    }

    @Override
    public double variance() {
        return Double.POSITIVE_INFINITY;
    }

    @Override
    public UnivariateFunction getProbabilityDensityFunction() {
        return null;
    }
}
