package dr.math.distributions;

import dr.math.UnivariateFunction;

/**
 * A distribution that is offset, so that the origin is greater than 0
 *
 * @author Alexei Drummond
 */
public class OffsetPositiveDistribution implements Distribution {

    /**
     * Constructor
     *
     * @param distribution distribution to offset
     * @param offset       a (positive) location parameter that allows this distribution to start
     *                     at a non-zero location
     */
    public OffsetPositiveDistribution(Distribution distribution, double offset) {

        if (offset < 0.0) throw new IllegalArgumentException();
        this.offset = offset;
        this.distribution = distribution;
    }

    /**
     * probability density function of the offset distribution
     *
     * @param x argument
     * @return pdf value
     */
    public final double pdf(double x) {
        if (offset < 0) return 0.0;
        return distribution.pdf(x - offset);
    }

    /**
     * log probability density function of the offset distribution
     *
     * @param x argument
     * @return pdf value
     */
    public final double logPdf(double x) {
        if (offset < 0) return Math.log(0.0);
        return distribution.logPdf(x - offset);
    }

    /**
     * cumulative density function of the offset distribution
     *
     * @param x argument
     * @return cdf value
     */
    public final double cdf(double x) {
        if (offset < 0) return 0.0;
        return distribution.cdf(x - offset);
    }

    /**
     * quantile (inverse cumulative density function) of the (offset) distribution
     *
     * @param y the p-value
     * @return icdf value
     */
    public final double quantile(double y) {
        return distribution.quantile(y) + offset;
    }

    /**
     * mean of the offset distribution
     *
     * @return mean
     */
    public final double mean() {
        return distribution.mean() + offset;
    }

    /**
     * variance of the offset distribution
     *
     * @return variance
     */
    public final double variance() {
        throw new UnsupportedOperationException();
    }

    public final UnivariateFunction getProbabilityDensityFunction() {
        return pdfFunction;
    }

    private UnivariateFunction pdfFunction = new UnivariateFunction() {
        public final double evaluate(double x) {
            return pdf(x);
        }

        public final double getLowerBound() {
            return offset + distribution.getProbabilityDensityFunction().getLowerBound();
        }

        public final double getUpperBound() {
            return offset + distribution.getProbabilityDensityFunction().getUpperBound();
        }
    };

    // the location parameter of the start of the positive distribution
    private double offset = 0.0;

    // the distribution to offset
    private Distribution distribution;
}
