package dr.geo.distributions;

import dr.geo.math.Space;

/**
 * @author Marc Suchard
 */
public class FisherBinghamKentDistribution extends HyperSphereDistribution {

    public FisherBinghamKentDistribution(int dim, Space space, double[] mean, double[] major, double[] minor,
                                         double kappa, double beta) {
        super(dim, space, mean, kappa);
        this.major = major;
        this.minor = minor;
        this.beta = beta;
    }

    public double logPdf(double[] x) {
        return logPdf(x, mean, major, minor);
    }

    public static double logPdf(double[] x, double[] gamma1, double[] gamma2, double[] gamma3) {
        return 0;
    }

    public String getType() {
        return "Fisher-Bingham (Kent)";
    }

    protected int getAllowableDim() {
        return 3;
    }

    private double[] major;
    private double[] minor;
    private double beta;
}