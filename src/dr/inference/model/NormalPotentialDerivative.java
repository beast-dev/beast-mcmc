package dr.inference.model;

/**
 * @author Max Tolkoff
 */
@Deprecated // TODO Should be implemented in NormalDistribution, etc.
public class NormalPotentialDerivative implements GradientWrtParameterProvider {
    double mean;
    double stdev;
    Parameter parameter;

    public NormalPotentialDerivative(double mean, double stdev, Parameter parameter){
        this.mean = mean;
        this.stdev = stdev;
        this.parameter = parameter;
    }

    @Override
    public Likelihood getLikelihood() {
        return null;
    }

    @Override
    public Parameter getParameter() {
        return parameter;
    }

    @Override
    public int getDimension() {
        return parameter.getDimension();
    }

    @Override
    public double[] getGradientLogDensity() {
        double[] derivative = new double[parameter.getDimension()];

        for (int i = 0; i < derivative.length; i++) {
            derivative[i] += (parameter.getParameterValue(i) - mean) / Math.sqrt(stdev);
            // TODO Should be?
            // derivative[i] = (mean - parameter.getParameterValue(i)) / variance;
        }

        return derivative;
    }
}
