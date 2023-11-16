package dr.inference.hmc;

import dr.evomodel.speciation.SpeciationLikelihoodGradient;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.math.MachineAccuracy;
import dr.math.matrixAlgebra.WrappedVector;
import dr.util.Transform;
import dr.xml.Reportable;

/**
 * @author Andy Magee
 * @author Yucai Shao
 */
public class GradientWrtIncrement implements GradientWrtParameterProvider, Reportable {

    private final GradientWrtParameterProvider gradient;
    private final Parameter incrementParameter;
    private final int dim;
    private final Transform incrementTransform;

    public GradientWrtIncrement(GradientWrtParameterProvider gradient, Parameter parameter, Transform incrementTransform) {
        this.gradient = gradient;
        this.incrementTransform = incrementTransform;
        this.incrementParameter = parameter;
        dim = gradient.getDimension();
    }

    @Override
    public Likelihood getLikelihood() {
        return gradient.getLikelihood();
    }

    @Override
    public Parameter getParameter() {
        return incrementParameter;
    }

    @Override
    public int getDimension() {
        return dim;
    }

    public double[] parameterFromIncrements(double[] delta) {
        double[] fx = new double[delta.length];
        fx[0] = delta[0];
        for (int i = 1; i < delta.length; i++) {
            fx[i] = fx[i-1] + delta[i];
        }
        return incrementTransform.inverse(fx, 0, delta.length);
    }

    @Override
    public double[] getGradientLogDensity() {
        // The gradient with respect to the variable-scale
        double[] grad = gradient.getGradientLogDensity();

        // The parameter on the scale of the gradient
        double[] gradScaleParameter = parameterFromIncrements(incrementParameter.getParameterValues());

        // The gradient with respect to the increments
        double[] incrementGrad = new double[dim];

        incrementGrad[dim - 1] = grad[dim - 1] * incrementTransform.gradient(gradScaleParameter[dim - 1]);
        for (int i = dim - 2; i > -1; i--) {
            incrementGrad[i] = grad[i] * incrementTransform.gradient(gradScaleParameter[i]) + incrementGrad[i + 1];
        }

        return incrementGrad;
    }

    @Override
    public String getReport() {
        StringBuilder sb = new StringBuilder();

        sb.append("Transform: ").append(incrementTransform.toString()).append("\n");
        sb.append("Gradient WRT increments: ").append(new dr.math.matrixAlgebra.Vector(getGradientLogDensity())).append("\n");
        sb.append("Gradient WRT parameters: ").append(new dr.math.matrixAlgebra.Vector(gradient.getGradientLogDensity())).append("\n");
        sb.append("Increments: ").append(new dr.math.matrixAlgebra.Vector(incrementParameter.getParameterValues())).append("\n");
        sb.append("Parameters: ").append(new dr.math.matrixAlgebra.Vector(parameterFromIncrements(incrementParameter.getParameterValues()))).append("\n");

        sb.append("Numerical gradient: ").append(new dr.math.matrixAlgebra.Vector(
                new CheckGradientNumerically(this, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, MachineAccuracy.SQRT_EPSILON, MachineAccuracy.SQRT_EPSILON).getNumericalGradient()
        )).append("\n");

        return sb.toString();
    }
}