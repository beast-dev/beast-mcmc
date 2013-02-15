package dr.math.distributions;

import dr.math.UnivariateFunction;
import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.AbstractContinuousDistribution;

/**
 * @author Andrew Rambaut
 */
public class TruncatedDistribution extends AbstractContinuousDistribution implements Distribution {

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

        if (source.getProbabilityDensityFunction().getLowerBound() > lower) {
            lower = source.getProbabilityDensityFunction().getLowerBound();
        }

        if (source.getProbabilityDensityFunction().getUpperBound() < upper) {
            upper = source.getProbabilityDensityFunction().getUpperBound();
        }

        this.lower = lower;
        this.upper = upper;

        if (!Double.isInfinite(this.lower)) {
            this.lowerCDF = source.cdf(lower);
        } else {
            this.lowerCDF = 0;
        }

        if (!Double.isInfinite(this.upper)) {
            this.normalization = source.cdf(upper) - lowerCDF;
        } else {
            this.normalization = 1.0 - lowerCDF;
        }
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
            cdf = 0.0;
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

        if (Double.isInfinite(lower) && Double.isInfinite(upper)) {
            return source.quantile(y);
        }

        try {
            return super.inverseCumulativeProbability(y);
        } catch (MathException e) {
//                throw MathRuntimeException.createIllegalArgumentException(                // AR - throwing exceptions deep in numerical code causes trouble. Catching runtime
            // exceptions is bad. Better to return NaN and let the calling code deal with it.
            return Double.NaN;

//                    "Couldn't calculate beta quantile for alpha = " + alpha + ", beta = " + beta + ": " +e.getMessage());
        }
    }

    /**
     * mean of the distribution
     *
     * @return mean
     */
    public double mean() {
        if (source != null) {
            return source.mean();
        } else {
            throw new IllegalArgumentException("Distribution is null");
        }
    }

    /**
     * variance of the distribution
     *
     * @return variance
     */
    public double variance() {
        if (source != null) {
            return source.variance();
        } else {
            throw new IllegalArgumentException("Distribution is null");
        }

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

    @Override
    protected double getInitialDomain(double v) {
        if (!Double.isInfinite(lower) && !Double.isInfinite(upper)) {
            return (upper + lower) / 2;
        } else if (!Double.isInfinite(upper)) {
            return upper / 2;
        } else if (!Double.isInfinite(lower)) {
            return lower * 2;
        }
        return v;
    }

    @Override
    protected double getDomainLowerBound(double v) {
        return lower;
    }

    @Override
    protected double getDomainUpperBound(double v) {
        return upper;
    }

    public double cumulativeProbability(double v) throws MathException {
        return cdf(v);
    }
}
