package dr.inference.hmc;

import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.math.matrixAlgebra.WrappedVector;
import dr.xml.Reportable;

/**
 * @author Andy Magee
 */
    public class NumericalGradient implements GradientWrtParameterProvider, Reportable {

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

    private final Likelihood likelihood;
    private final Parameter parameter;

    @Override
    public String getReport() {
        long start = System.currentTimeMillis();
        double[] result = getGradientLogDensity();
        long duration = System.currentTimeMillis() - start;

        return "Numeric gradient = " + new WrappedVector.Raw(result) + " (in " + duration + "ms)";
    }
}
