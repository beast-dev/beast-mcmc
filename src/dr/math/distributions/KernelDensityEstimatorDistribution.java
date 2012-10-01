package dr.math.distributions;

import dr.math.UnivariateFunction;

/**
 * @author Marc Suchard
 */
public abstract class KernelDensityEstimatorDistribution implements Distribution {

    public KernelDensityEstimatorDistribution(Double[] sample, Double lowerBound, Double upperBound, Double bandWidth) {

        this.sample = new double[sample.length];
        for (int i = 0; i < sample.length; i++) {
            this.sample[i] = sample[i];
        }
        this.N = sample.length;
        processBounds(lowerBound, upperBound);
        setBandWidth(bandWidth);
    }

    abstract protected double evaluateKernel(double x);

    abstract protected void processBounds(Double lowerBound, Double upperBound);

    abstract protected void setBandWidth(Double bandWidth);

    /**
     * probability density function of the distribution
     *
     * @param x argument
     * @return pdf value
     */
    public double pdf(double x) {
        return evaluateKernel(x);
    }

    /**
     * the natural log of the probability density function of the distribution
     *
     * @param x argument
     * @return log pdf value
     */
    public double logPdf(double x) {
        return Math.log(pdf(x));
    }

    /**
     * cumulative density function of the distribution
     *
     * @param x argument
     * @return cdf value
     */
    public double cdf(double x) {
        throw new RuntimeException("Not Implemented.");
    }

    /**
     * quantile (inverse cumulative density function) of the distribution
     *
     * @param y argument
     * @return icdf value
     */
    public double quantile(double y) {
        throw new RuntimeException("Not Implemented.");
    }

    /**
     * mean of the distribution
     *
     * @return mean
     */
    public double mean() {
        throw new RuntimeException("Not Implemented.");
    }

    /**
     * variance of the distribution
     *
     * @return variance
     */
    public double variance() {
        throw new RuntimeException("Not Implemented.");
    }

    /**
     * @return a probability density function representing this distribution
     */
    public UnivariateFunction getProbabilityDensityFunction() {
        throw new RuntimeException("Not Implemented.");
    }

    public double getBandWidth() {
        return bandWidth;
    }

    public enum Type {
        GAUSSIAN("Gaussian"),
        GAMMA("Gamma"),
        LOGTRANSFORMEDGAUSSIAN("LogTransformedGaussian"),
        BETA("Beta");

        private Type(String text) {
            this.text = text;
        }

        public String getText() {
            return text;
        }

        public static Type parseFromString(String text) {
            for (Type format : Type.values()) {
                if (format.getText().compareToIgnoreCase(text) == 0)
                    return format;
            }
            return null;
        }

        private final String text;
    }

    protected int N;
    protected double lowerBound;
    protected double upperBound;
    protected double bandWidth;
    protected double[] sample;
}
