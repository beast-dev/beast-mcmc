package dr.inference.hmc;

import dr.evomodel.speciation.SpeciationLikelihoodGradient;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.math.MachineAccuracy;
import dr.math.matrixAlgebra.WrappedVector;
import dr.xml.Reportable;

/**
 * @author Andy Magee
 */
public class GradientWrtIncrement implements GradientWrtParameterProvider, Reportable {

    private final GradientWrtParameterProvider gradient;
    private final Parameter incrementParameter;
    private final int dim;
    private final IncrementTransformType type;

    public GradientWrtIncrement(GradientWrtParameterProvider gradient, Parameter parameter, IncrementTransformType type) {
        this.gradient = gradient;
        this.type = type;
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

    private double getGradientOfTransformation(double y) {
        if (this.type.getTransformType().equalsIgnoreCase("log")) {
            return y;
        } else if (this.type.getTransformType().equalsIgnoreCase("logit")) {
            LogitTransform logitType = (LogitTransform) this.type;
            double scaledY = (y - logitType.getLower()) / (logitType.getUpper() - logitType.getLower());
            return (logitType.getUpper() - logitType.getLower()) * scaledY * (1-scaledY);
        } else {
            throw new RuntimeException("Not implemented!");
        }
    }
    @Override
    public double[] getGradientLogDensity() {
        // The gradient with respect to the variable-scale
        double[] grad = gradient.getGradientLogDensity();

        // The parameter on the scale of the gradient
        double[] gradScaleParameter = type.parameterFromIncrements(incrementParameter.getParameterValues());

        // The gradient with respect to the increments
        double[] incrementGrad = new double[dim];

        // TODO: is this only right for log-transforms?
        incrementGrad[dim - 1] = grad[dim - 1] * getGradientOfTransformation(gradScaleParameter[dim - 1]);
        for (int i = dim - 2; i > -1; i--) {
            incrementGrad[i] = grad[i] * getGradientOfTransformation(gradScaleParameter[i]) + incrementGrad[i + 1];
        }

        return incrementGrad;
    }

    @Override
    public String getReport() {
        StringBuilder sb = new StringBuilder();

        sb.append("Gradient WRT increments: ").append(new dr.math.matrixAlgebra.Vector(getGradientLogDensity())).append("\n");
        sb.append("Gradient WRT parameters: ").append(new dr.math.matrixAlgebra.Vector(gradient.getGradientLogDensity())).append("\n");
        sb.append("Increments: ").append(new dr.math.matrixAlgebra.Vector(incrementParameter.getParameterValues())).append("\n");
        sb.append("Parameters: ").append(new dr.math.matrixAlgebra.Vector(type.parameterFromIncrements(incrementParameter.getParameterValues()))).append("\n");

        sb.append("Numerical gradient: ").append(new dr.math.matrixAlgebra.Vector(
                new CheckGradientNumerically(this, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, MachineAccuracy.SQRT_EPSILON, MachineAccuracy.SQRT_EPSILON).getNumericalGradient()
        )).append("\n");

        return sb.toString();
    }
}
