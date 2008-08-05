package dr.math;

/**
 * @author Alexei Drummond
 */
public class LaplaceDistribution implements Distribution {

    //Parameters
    double mu;

    //Parameters
    double beta;

    //Parameters
    double c;

    /**
     * This general constructor creates a new laplace distribution with a
     * specified rate
     */
    public LaplaceDistribution(double mu, double beta) {
        setParameters(mu, beta);
    }

    /**
     * This default constructor creates a new laplace distribution with rate 1
     */
    public LaplaceDistribution() {
        this(0, 1);
    }

    /**
     * Set parameters and assign the default partition
     */
    public void setParameters(double k, double b) {
        if (b <= 0) b = 1;
        mu = k;
        beta = b;

        //Normalizing constant
        c = 1 / (2 * beta);
    }

    /**
     * Get center parameters
     */
    public double getMu() {
        return mu;
    }

    /**
     * Get scale parameters
     */
    public double getBeta() {
        return beta;
    }

    /**
     * Maximum value of getDensity function
     */
    public double getMaxDensity() {
        return c;
    }

    /**
     * Cumulative distribution function
     */
    public double cdf(double x) {
        if (x == mu) return 0.5;
        else return (0.5) * (1 + ((x - mu) / Math.abs(x - mu))
                * (1 - Math.exp(-Math.abs(x - mu) / beta)));
    }

    /**
     * Density function
     */
    public double pdf(double x) {
        return c * Math.exp(-Math.abs(x - mu) / beta);
    }

    public double logPdf(double x) {
        return Math.log(c) - (Math.abs(x - mu) / beta);
    }

    public double quantile(double y) {
        throw new UnsupportedOperationException();
    }

    public double mean() {
        return mu;
    }

    public double variance() {
        return 2 * beta * beta;
    }

    public UnivariateFunction getProbabilityDensityFunction() {
        throw new UnsupportedOperationException();
    }
}
