package dr.evomodel.continuous;

import dr.inference.model.Parameter;

/**
 * @author Marc Suchard
 *
 *
 * Follows a bifractional diffusion model developed in:
 * Brockman D, Hufnagel L and Geisel T (2006) The scaling laws of human travel. Nature 439, 462 - 465
 *
 *
 */
public class BifractionalDiffusionModel extends MultivariateDiffusionModel {

    public BifractionalDiffusionModel(Parameter alpha, Parameter beta) {
        super();
        this.alpha = alpha;
        this.beta = beta;
        addVariable(alpha);
        addVariable(beta);
    }


    protected double calculateLogDensity(double[] start, double[] stop, double time) {
        // Compute finite-time transition probability

        // Equation (2) from Brockman, Hufnagel and Geisel (2006)
        final double ratio = alpha.getParameterValue(0) / beta.getParameterValue(0);
        final double r = distanceEuclidean(start, stop);
        final double scaledTime = Math.pow(time,ratio);        
        return -ratio * Math.log(time) + logUniversalScalingFunction(r / scaledTime);
    }

    private double distanceEuclidean(double[] start, double[] stop) {
        final int dim = start.length;
        double total = 0;
        for(int i=0; i<dim; i++) {
            final double dX = stop[i] - start[i];
            total += dX*dX;
        }
        return Math.sqrt(total);
    }

    private double logUniversalScalingFunction(double x) {
        return x - x;
    }

    protected void calculatePrecisionInfo() {
        // Precompute normalizing constants if necessary
    }

    private Parameter alpha;
    private Parameter beta;

}
