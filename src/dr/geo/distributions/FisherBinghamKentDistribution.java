package dr.geo.distributions;

/**
 * @author Marc Suchard
 */
public class FisherBinghamKentDistribution extends HyperSphereDistribution {

    public FisherBinghamKentDistribution(int dim, double[] mean, double[] major, double[] minor) {
        super(dim, mean);
        this.major = major;
        this.minor = minor;
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
}