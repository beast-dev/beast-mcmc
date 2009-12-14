package dr.geo.distributions;

/**
 * @author Marc Suchard
 */
public class VonMisesDistribution extends HyperSphereDistribution {

    public VonMisesDistribution(int dim, double[] mean) {
        super(dim, mean);
    }

    public double logPdf(double[] x) {
        return 0;
    }

    public String getType() {
        return "von Mises";
    }

    protected int getAllowableDim() {
        return 2;
    }
}