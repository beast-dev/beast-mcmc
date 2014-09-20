package dr.math.distributions;

import dr.math.ErrorFunction;
import dr.math.UnivariateFunction;

/**
 * User: Michael Defoin Platel
 * Date: Oct 31, 2007
 * Time: 2:32:37 PM
 */
public class TruncatedNormalDistribution implements Distribution {

    public double getMean() {
        return m;
    }

    public void setMdean(double m) {
        this.m = m;
    }

    public double getSD() {
        return sd;
    }

    public void setSD(double sd) {
        this.sd = sd;
    }

    public double getLower() {
        return lower;
    }

    public void setLower(double lower) {
        this.lower = lower;
    }

    public double getUpper() {
        return upper;
    }

    public void setUpper(double upper) {
        this.upper = upper;
    }

    public TruncatedNormalDistribution(double mean, double sd, double lower, double upper) {


        if (lower == upper)
            upper += 1.E-4;

        if (sd == 0.)
            sd = 1.E-5;

        this.m = mean;
        this.sd = sd;
        this.lower = lower;
        this.upper = upper;

        double upperCDF;
        double lowerCDF;

        if(upper!=Double.POSITIVE_INFINITY){
            upperCDF = standardNormalCdf((upper-mean)/sd);
        } else {
            upperCDF = 1;
        }

        if(lower!=Double.NEGATIVE_INFINITY){
            lowerCDF = standardNormalCdf((lower-mean)/sd);
        } else {
            lowerCDF = 0;
        }

        this.T = upperCDF - lowerCDF;
    }


    public double pdf(double x) {
        if (x >= upper || x < lower)
            return 0.0;
        else
            return (standardNormalPdf((x - m) / sd) / sd) / T;
    }

    public double logPdf(double x) {
        return Math.log(pdf(x));
    }

    public double cdf(double x) {
        double cdf;
        if (x < lower)
            cdf = 0.;
        else if (x >= lower && x < upper)
            if(lower!=Double.NEGATIVE_INFINITY){
                cdf = (standardNormalCdf((x - m) / sd) - standardNormalCdf((lower - m) / sd)) / T;
            } else {
                cdf = (standardNormalCdf((x - m) / sd)) / T;
            }
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

    public double mean() {
        return mean(m, sd, lower, upper);
    }

    public double variance() {
        return mean(m, sd, lower, upper);
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

    /**
     * probability density function  of the standard normal distribution
     *
     * @param x argument
     * @return pdf at x
     */
    public static double standardNormalPdf(double x) {
        double a = 1.0 / (Math.sqrt(2.0 * Math.PI));
        double b = -(x) * (x) / (2.0);

        return a * Math.exp(b);
    }

    /**
     * the natural log of the probability density function of the standard normal distribution
     *
     * @param x argument
     * @return log pdf at x
     */
    public static double logStandardNormalPdf(double x) {
        /*Can throw an expcetion if x not in the range [lower,upper]*/
        return Math.log(standardNormalPdf(x));
    }

    /**
     * cumulative density function of the standard normal distribution
     *
     * @param x argument
     * @return cdf at x
     */
    public static double standardNormalCdf(double x) {
        double a = (x) / (Math.sqrt(2.0));

        return 0.5 * (1.0 + ErrorFunction.erf(a));
    }

    /**
     * mean
     *
     * @param m     mean
     * @param sd    standard deviation
     * @param lower the lower limit
     * @param upper the upper limit
     * @return mean
     */
    public static double mean(double m, double sd, double lower, double upper) {
        double au, al, pu, pl, cu, cl;

        au = (upper - m) / sd;
        al = (lower - m) / sd;

        pu = standardNormalPdf(au);
        pl = standardNormalPdf(al);

        cu = standardNormalCdf(au);
        cl = standardNormalCdf(al);

        return m - sd * (pu - pl) / (cu - cl);
    }

    /**
     * variance
     *
     * @param m     mean
     * @param sd    standard deviation
     * @param lower the lower limit
     * @param upper the upper limit
     * @return variance
     */
    public static double variance(double m, double sd, double lower, double upper) {
        double au, al, pu, pl, cu, cl, T1, T2;

        au = (upper - m) / sd;
        al = (lower - m) / sd;

        pu = standardNormalPdf(au);
        pl = standardNormalPdf(al);

        cu = standardNormalCdf(au);
        cl = standardNormalCdf(al);

        T1 = (au * pu - al * pl) / (cu - cl);
        T2 = (pu - pl) / (cu - cl);

        return sd * sd * (1.0 - T1 - T2 * T2);
    }


    private double m;
    private double sd;
    private double lower;
    private double upper;
    private double T;
}
