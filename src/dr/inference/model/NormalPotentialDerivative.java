package dr.inference.model;

/**
 * @author Max Tolkoff
 */
public class NormalPotentialDerivative implements GradientProvider {
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
        throw new RuntimeException("Not yet implemented");
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
        }

        return derivative;
    }
}
