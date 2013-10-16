package dr.math.distributions;

/**
 * @author Marc A. Suchard
 */
public interface MultivariateDistribution {

    public double logPdf(double[] x);

    public double[][] getScaleMatrix();

    public double[] getMean();

    public String getType();

}
