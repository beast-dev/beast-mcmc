package dr.math.distributions;

import dr.math.ErrorFunction;
import dr.math.UnivariateFunction;

/**
 * @author Andrew Rambaut
 */
public class TruncatedDistribution implements Distribution {

    public double getLower() {
        return lower;
    }

    public double getUpper() {
        return upper;
    }

    public TruncatedDistribution(Distribution source, double lower, double upper) {
        this.source = source;

        if (lower == upper) {
            throw new IllegalArgumentException("upper equals lower");
        }

        this.lower = lower;
        this.upper = upper;

        this.lowerCDF = source.cdf(lower);
        this.normalization = source.cdf(upper) - lowerCDF;
    }


    public double pdf(double x) {
        if (x >= upper && x < lower)
            return 0.0;
        else
            return source.pdf(x) / normalization;
    }

    public double logPdf(double x) {
        return Math.log(pdf(x));
    }

    public double cdf(double x) {
        double cdf;
        if (x < lower)
            cdf = 0.;
        else if (x >= lower && x < upper)
            cdf = (source.cdf(x) - lowerCDF) / normalization;
        else
            cdf = 1.0;

        return cdf;
    }

    public double quantile(double y) {

        if (y == 0)
            return lower;

        if (y == 1.0)
            return upper;

        return quantileSearch(y, lower, upper, 20);
    }

    /*Implements a geometic search for the quantiles*/
    private double quantileSearch(double y, double l, double u, int step) {
        double q, a;

        q = (u + l) / 2.0;

        if (step == 0 || q == l || q == u)
            return q;

        a = cdf(q);

        if (y <= a)
            return quantileSearch(y, l, q, step - 1);
        else
            return quantileSearch(y, q, u, step - 1);
    }

    /**
     * mean of the distribution
     *
     * @return mean
     */
    public double mean() {
        throw new UnsupportedOperationException("Not Implemented.");
    }

    /**
     * variance of the distribution
     *
     * @return variance
     */
    public double variance() {
        throw new UnsupportedOperationException("Not Implemented.");
    }

    public UnivariateFunction getProbabilityDensityFunction() {
        return pdfFunction;
    }

    private UnivariateFunction pdfFunction = new UnivariateFunction() {
        public final double evaluate(double x) {
            return pdf(x);
        }

        public final double getLowerBound() {
            return lower;
        }

        public final double getUpperBound() {
            return upper;
        }
    };

    private final Distribution source;
    private final double lower;
    private final double upper;
    private final double normalization;
    private final double lowerCDF;
}
