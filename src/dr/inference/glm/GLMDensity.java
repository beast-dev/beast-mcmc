package dr.inference.glm;

import dr.math.distributions.NormalDistribution;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public interface GLMDensity {

    double getLogPDF(double[] beta, double[] X);

}
