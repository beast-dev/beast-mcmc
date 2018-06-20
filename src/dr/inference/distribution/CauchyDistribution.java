package dr.inference.distribution;



import dr.math.UnivariateFunction;
import dr.math.distributions.Distribution;
import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.AbstractContinuousDistribution;
import org.apache.commons.math.distribution.CauchyDistributionImpl;

public class CauchyDistribution extends AbstractContinuousDistribution implements Distribution {
    CauchyDistributionImpl distribution;
    double median;
    double scale;

    public CauchyDistribution(double median, double scale){
        distribution = new CauchyDistributionImpl(median, scale);
        this.median = median;
        this.scale = scale;
    }


    @Override
    public double pdf(double x) {
        return distribution.density(x);
    }

    @Override
    public double logPdf(double x) {
        return Math.log(distribution.density(x));
    }

    @Override
    public double cdf(double x) {
        return distribution.cumulativeProbability(x);
    }

    @Override
    public double quantile(double y) {
        return distribution.inverseCumulativeProbability(y);
    }

    @Override
    public double mean() {
        return Double.NaN;
    }

    @Override
    public double variance() {
        return Double.POSITIVE_INFINITY;
    }

    @Override
    public final UnivariateFunction getProbabilityDensityFunction() {
        return pdfFunction;
    }

    private final UnivariateFunction pdfFunction = new UnivariateFunction() {
        public final double evaluate(double x) {
            return pdf(x);
        }

        public final double getLowerBound() {
            return Double.NEGATIVE_INFINITY;
        }

        public final double getUpperBound() {
            return Double.POSITIVE_INFINITY;
        }
    };

    @Override
    protected double getInitialDomain(double v) {
        return Double.POSITIVE_INFINITY;
    }

    @Override
    protected double getDomainLowerBound(double v) {
        return Double.NEGATIVE_INFINITY;
    }

    @Override
    protected double getDomainUpperBound(double v) {
        return Double.POSITIVE_INFINITY;
    }

    @Override
    public double cumulativeProbability(double v) throws MathException {
        return 0;
    }
}
