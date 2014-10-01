package dr.geo.distributions;

import cern.jet.math.Bessel;
import dr.geo.math.Space;

/**
 * @author Marc Suchard
 */
public class VonMisesDistribution extends HyperSphereDistribution {
        
    public VonMisesDistribution(int dim, Space space, double[] mean, double kappa) {
        super(dim, space, mean, kappa);
    }

    public double logPdf(double[] x) {
        return logPdf(x, mean, kappa, space);
    }

    public String getType() {
        return "von Mises";
    }

    protected int getAllowableDim() {
        return 2;
    }

    private static double logNormalizationConstant(double kappa) {
         return -LOG_2_PI - Math.log(Bessel.i0(kappa));
    }

    public static double logPdf(double[] x, double[] mean, double kappa, Space space) {
        return logNormalizationConstant(kappa) + kappa * HyperSphereDistribution.innerProduct(x, mean, space);
    }

    public static void main(String[] arg) {        
        // Test in radians
        double kappa = 2;
        double[] mean = { 1 };
        double[] x = { -4 };
        System.err.println("logP = "+logPdf(x, mean, kappa, Space.RADIANS)+" ?= -2.094546");
    }
}