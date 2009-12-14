package dr.geo.distributions;

import dr.math.distributions.MultivariateDistribution;

/**
 * @author Marc Suchard
 *
 * Abstract base class for probability distributions on a hyper-sphere S^{d-1} in R^{d}.
 * Usual cases are: circle (d=2) and sphere (d=3)
 */

public abstract class HyperSphereDistribution implements MultivariateDistribution {

    public HyperSphereDistribution(int dim, double[] mean) {
        this.dim = dim;
        this.mean = mean;
        checkImplementation();
    }

    public abstract double logPdf(double[] x);

    public double[][] getScaleMatrix() {
        return new double[0][];
    }

    public double[] getMean() {
        return new double[0];
    }

    public String getType() {
        return "Hyperspherical";
    }

    protected int getAllowableDim() {
        return dim;
    }

    protected void checkImplementation() {
        final int allowableDim = getAllowableDim();
        if (allowableDim != dim) {
            throw new RuntimeException(getType() + " distribution is not implemented in R^"+dim);
        }
    }

    protected int dim;
    protected double[] mean;
}
