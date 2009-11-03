package dr.math.distributions;

import dr.math.Poisson;
import dr.math.UnivariateFunction;
import org.apache.commons.math.MathException;

/**
 * @author Alexei Drummond
 * @version $Id$
 */
public class PoissonDistribution implements Distribution {

    org.apache.commons.math.distribution.PoissonDistribution distribution;

    public PoissonDistribution(double mean) {
        distribution = new org.apache.commons.math.distribution.PoissonDistributionImpl(mean);
    }

    public double pdf(double x) {
        return distribution.probability(x);
    }

    public double logPdf(double x) {

        double pdf = distribution.probability(x);
        if (pdf == 0 || Double.isNaN(pdf)) { // bad estimate
            final double mean = mean();
            return x * Math.log(mean) - Poisson.gammln(x + 1) - mean;
        }
        return Math.log(pdf);

    }

    public double cdf(double x) {
        try {
            return distribution.cumulativeProbability(x);
        } catch (MathException e) {
            throw new RuntimeException(e);
        }
    }

    public double quantile(double y) {
        try {
            return distribution.inverseCumulativeProbability(y);
        } catch (MathException e) {
            throw new RuntimeException(e);
        }
    }

    public double mean() {
        return distribution.getMean();
    }

    public double variance() {
        return distribution.getMean();
    }

    public UnivariateFunction getProbabilityDensityFunction() {
        throw new RuntimeException();
    }

    public double truncatedMean(int max) {

        double CDF = 0;
        double mean = 0;
        for(int i=0; i<=max; i++) {
            double p = distribution.probability(i);
            mean += i*p;
            CDF += p;
        }
        return mean / CDF;        
    }
}
