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
        return new CheckGradientNumerically(this,
                Double.NEGATIVE_INFINITY,
                Double.POSITIVE_INFINITY,
                0.0,
                0.0).getNumericalGradient();
    }

    private Likelihood likelihood = null;
    private Parameter parameter = null;
}
