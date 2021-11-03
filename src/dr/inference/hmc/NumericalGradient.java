package dr.inference.hmc;

import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.math.MachineAccuracy;

/**
 * @author Andy Magee
 */
public class NumericalGradient implements GradientWrtParameterProvider {

    public NumericalGradient(Likelihood likelihood, Parameter parameter) {
        this.parameter = parameter;
        this.likelihood = likelihood;
    }

    @Override
    public Likelihood getLikelihood() {
        return likelihood;
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
        final int dim = parameter.getDimension();
        double[] gradient = new double[dim];

        final double[] oldValues = parameter.getParameterValues();

        double[] likelihoodPlus = new double[dim];
        double[] likelihoodMinus = new double[dim];

        double[] h = new double[dim];
        for (int i = 0; i < dim; i++) {
            h[i] = MachineAccuracy.SQRT_SQRT_EPSILON * (Math.abs(oldValues[i]) + 1.0);
            parameter.setParameterValue(i, oldValues[i] + h[i]);
            likelihoodPlus[i] = likelihood.getLogLikelihood();

            parameter.setParameterValue(i, oldValues[i] - h[i]);
            likelihoodMinus[i] = likelihood.getLogLikelihood();
            parameter.setParameterValue(i, oldValues[i]);
        }

        for (int i = 0; i < dim; i++) {
            gradient[i] = (likelihoodPlus[i] - likelihoodMinus[i]) / (2.0 * h[i]);
        }

        return gradient;
    }

    private Likelihood likelihood = null;
    private Parameter parameter = null;
}
