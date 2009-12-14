package dr.geo.distributions;

/**
 * @author Marc Suchard
 */
public class VonMisesFisherDistribution extends HyperSphereDistribution {

    public VonMisesFisherDistribution(int dim, double[] mean) {
        super(dim, mean);
    }

    public double logPdf(double[] x) {
        return 0;
    }

    public String getType() {
        return "von Mises-Fisher";
    }

    protected int getAllowableDim() {
        return 3;
    }
}
