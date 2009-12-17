package dr.math.distributions;

import dr.math.MittagLefflerFunction;
import dr.math.UnivariateFunction;

/**
 * @author Marc Suchard
 */

public class AsymptoticBifractionalDiffusionDensity implements Distribution {

    public AsymptoticBifractionalDiffusionDensity() {


    }

    


    /**
     * probability density function of the distribution
     *
     * @param x argument
     * @return pdf value
     */
    public double pdf(double x) {
        throw new RuntimeException("Not yet implemented");
    }

    /**
     * the natural log of the probability density function of the distribution
     *
     * @param x argument
     * @return log pdf value
     */
    public double logPdf(double x) {
        throw new RuntimeException("Not yet implemented");
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


    private double alpha;
    private double beta;
    private MittagLefflerFunction mlFunc;
}
