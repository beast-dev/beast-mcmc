package dr.math.distributions;

import dr.math.UnivariateFunction;

/**
 * Created by mandevgill on 11/15/14.
 */
public class DiscreteUniformDistribution implements Distribution {

    public DiscreteUniformDistribution(double lower, double upper) {
        this.lower = lower;
        this.upper = upper;
        this.n = upper - lower + 1;
    }

    public double pdf(double x) {
        if (x < lower) {
            return 0;
        } else if (x > upper) {
            return 0;
        } else {
            return 1 / n;
        }
    }

    public double logPdf(double x) {
        /*
        if (x < 0)  return Double.NEGATIVE_INFINITY;
        double r = -1 * (mean*mean) / (mean - stdev*stdev);
        double p = mean / (stdev*stdev);
        return Math.log(Math.pow(1-p,x)) + Math.log(Math.pow(p, r)) + GammaFunction.lnGamma(r + x) - GammaFunction.lnGamma(r) - GammaFunction.lnGamma(x+1);
        */
        if (x < lower) {
            return Double.NEGATIVE_INFINITY;
        } else if (x > upper) {
            return Double.NEGATIVE_INFINITY;
        } else {
            return -Math.log(n);
        }
    }

    public double cdf(double x) {
        return (x - lower + 1) / n;
    }

    public double quantile(double y) {
        // fill in
        return Double.NaN;
    }

    public double mean() {
        return (lower + upper) / 2;
    }

    public double variance() {
        return (n * n - 1) / 12;
    }

    public UnivariateFunction getProbabilityDensityFunction() {
        throw new RuntimeException();
    }

    double lower;
    double upper;
    double n;

}
