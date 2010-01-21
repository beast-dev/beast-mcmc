package dr.math.distributions;

import dr.math.UnivariateFunction;
import cern.jet.stat.Gamma;

/**
 * @author Marc Suchard
 */

public class AsymptoticBifractionalDiffusionDensity implements Distribution {

    public AsymptoticBifractionalDiffusionDensity(double alpha, double beta) {
        this.alpha = alpha;
        this.beta = beta;
    }

    /**
     * probability density function of the distribution
     *
     * @param x argument
     * @return pdf value
     */
    public double pdf(double x) {
        return Math.exp(logPdf(x));
    }

    /**
     * the natural log of the probability density function of the distribution
     *
     * @param x argument
     * @return log pdf value
     */
    public double logPdf(double x) {
        return logPdf(x, 1.0, alpha, beta);
    }

    /**
     * cumulative density function of the distribution
     *
     * @param x argument
     * @return cdf value
     */
    public double cdf(double x) {
        throw new RuntimeException("Not yet implemented");
    }

    /**
     * quantile (inverse cumulative density function) of the distribution
     *
     * @param y argument
     * @return icdf value
     */
    public double quantile(double y) {
        throw new RuntimeException("Not yet implemented");
    }

    /**
     * mean of the distribution
     *
     * @return mean
     */
    public double mean() {
        throw new RuntimeException("Not yet implemented");
    }

    /**
     * variance of the distribution
     *
     * @return variance
     */
    public double variance() {
        throw new RuntimeException("Not yet implemented");
    }

    /**
     * @return a probability density function representing this distribution
     */
    public UnivariateFunction getProbabilityDensityFunction() {
        throw new RuntimeException("Not yet implemented");
    }

   /**
     * Taken from:  Saichev AI and Zaslavsky GM (1997) Fractional kinetic equations: solutions and applications.
     *              Chaos, 7, 753-764
     * @param x evaluation point
     * @param t evaluation time
     * @param alpha coefficient
     * @param beta coefficient
     * @return probability density
     */

    public static double logPdf(double x, double t, double alpha, double beta) {
        final double mu = beta / alpha;
        final double absX = Math.abs(x);
        final double absY = absX / Math.pow(t,mu);
        double density = 0;
        double incr = Double.MAX_VALUE;
        int m = 1; // 0th term = 0 \propto cos(pi/2)
        int sign = -1;

        while (incr > eps && m < max) {
            incr =  sign / Math.pow(absY, m * alpha)
                         * Gamma.gamma(m * alpha + 1)
                         / Gamma.gamma(m * beta + 1)
                         * Math.cos(halfPI * (m * alpha + 1));
            density += incr;
            sign *= -1;
            m++;
        }

        return Math.log(density / (Math.PI * absX));
    }

    private static double eps = 1E-10;
    private static int max = 1000;
    private static double halfPI = Math.PI / 2.0;

    private double alpha;
    private double beta;
}
