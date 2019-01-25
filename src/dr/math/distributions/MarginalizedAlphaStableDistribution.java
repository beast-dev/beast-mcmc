package dr.math.distributions;

import dr.math.UnivariateFunction;

/**
 * @author Marc A. Suchard
 * @author Philippe Lemey
 */
public class MarginalizedAlphaStableDistribution implements Distribution {

    private final double scale;
    private final double alpha;

    public MarginalizedAlphaStableDistribution(double scale, double alpha) {
        this.scale = scale;
        this.alpha = alpha;
    }

    @Override
    public double pdf(double x) {
        return Math.exp(logPdf(x));
    }

    @Override
    public double logPdf(double x) {
        return logPdf(x, scale, alpha);
    }

    @Override
    public double cdf(double x) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double quantile(double y) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double mean() {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double variance() {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public UnivariateFunction getProbabilityDensityFunction() {
        throw new RuntimeException("Not yet implemented");
    }

    public static double logPdf(double x, double scale, double alpha) {
        return -Math.log(scale) - alpha * Math.abs(x) / scale;
    }

    public static double gradLogPdf(double x, double scale, double alpha) {
        throw new RuntimeException("Not yet implemented");
    }
}
