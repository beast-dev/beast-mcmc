package dr.inference.model;

import dr.inference.distribution.DistributionLikelihood;
import dr.math.distributions.NormalDistribution;

/**
 * @author Max Tolkoff
 */
public class NormalPotentialDerivative implements PotentialDerivativeInterface {
    double mean;
    double stdev;
    Parameter parameter;

    public NormalPotentialDerivative(double mean, double stdev, Parameter parameter){
        this.mean = mean;
        this.stdev = stdev;
        this.parameter = parameter;
    }

    @Override
    public double[] getDerivative() {
        double[] derivative = new double[parameter.getDimension()];

        for (int i = 0; i < derivative.length; i++) {
            derivative[i] += (parameter.getParameterValue(i) - mean) / Math.sqrt(stdev);
        }

        return derivative;
    }
}
