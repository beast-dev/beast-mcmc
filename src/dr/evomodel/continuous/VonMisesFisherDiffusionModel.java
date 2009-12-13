package dr.evomodel.continuous;

import dr.inference.model.Parameter;
import dr.evolution.continuous.SphericalPolarCoordinates;
import cern.jet.math.Bessel;

/**
 * @author Marc Suchard
 *         <p/>
 *         The Von Mises-Fisher distributions provides an approximation to diffusion on a sphere
 *         <p/>
 *         Main paper:
 *         R.A. Fisher (1953) Dispersion on a sphere. Proceedings of the Royal Society of London, Series A, 217, 295 - 305.
 *         <p/>
 *         Generalizations are discussed in:
 *         K. Shimizu and K. Iida, Pearson type VII distributions on spheres, Commun. Stat.ÐTheory Meth. 31 (2002), pp. 513Ð526
 *         K. shimizu and H.Y. Slew, The generalized symmetric Laplace distribution on the sphere, Statistical Methology, 5, 2008
 */

public class VonMisesFisherDiffusionModel extends MultivariateDiffusionModel {

    private static double LOG_2_PI = Math.log(2 * Math.PI);


    public VonMisesFisherDiffusionModel(Parameter concentrationParameter) {
        this(3, concentrationParameter); // Default is distribution on a sphere
    }

    public VonMisesFisherDiffusionModel(int p, Parameter concentrationParameter) {

        super();
        this.p = p;
        this.concentrationParameter = concentrationParameter;
        calculatePrecisionInfo();
        addVariable(concentrationParameter);
    }

    protected void calculatePrecisionInfo() {
        if (p == 1 || p > 3) {
            throw new RuntimeException("Von Mises-Fisher distribution only implemented for circles and spheres");
        }
        // Nothing gets stored
    }

    private static double computeLogNormalization(double kappa, int p) {

        if (p == 3) {
            // 'sinh' has some numerical instability for small arguments
            if (kappa < 1E-10) {
                return -Math.log(2) - LOG_2_PI;
            }
            return Math.log(kappa) - LOG_2_PI - Math.log(Math.exp(+kappa) - Math.exp(-kappa));
        } else if (p == 2) { // Bessel function of order (p/2-1)
            return Bessel.i0(kappa);
        }
        return 0;
    }

    private double cartesianInnerProduct(double[] x, double[] y, double radius) {

        // x[] and y[] should be in the form (lat, long)
        if (x.length != 2 || y.length != 2) {
            throw new RuntimeException("Wrong dimensions");
        }

        final SphericalPolarCoordinates coordX = new SphericalPolarCoordinates(x[0], x[1], radius);
        final SphericalPolarCoordinates coordY = new SphericalPolarCoordinates(y[0], y[1], radius);

        return coordX.getCartesianCoordinates().dot(coordY.getCartesianCoordinates());
    }

    protected double calculateLogDensity(double[] start, double[] stop, double time) {

        double innerProduct;
        if (p == 3 && start.length == 2) { // Given in (lat, long)
            innerProduct = cartesianInnerProduct(start, stop, 1.0);   // Distributional form assume ||x|| = ||y|| = 1
        } else if (p == 2 && start.length == 1) {
            innerProduct = Math.cos(stop[0] - start[0]); // Assumed to already be in radians   
        } else {
            innerProduct = 0;
            for(int i = 0; i < start.length; i++) {
                innerProduct += start[i] * stop[i];
            }
        }
        final double kappa = concentrationParameter.getParameterValue(0) / time;
        return computeLogNormalization(kappa,p) + kappa * innerProduct;
    }

    private Parameter concentrationParameter;
    private int p;

    public static void main(String[] arg) {

        double kappa = 1E-15;
        double r = computeLogNormalization(kappa,3);
        System.err.println("t(r) = " + r);
    }

}
