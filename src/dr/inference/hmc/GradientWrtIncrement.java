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

    public enum IncrementTransformType {
        LOG("log") {
            public double getDerivativeOfInverseTransform(double x) {
                return Math.exp(x);
            }

            public double[] parameterFromIncrements(double[] delta) {
                double[] fx = new double[delta.length];
                fx[0] = delta[0];
                for (int i = 1; i < delta.length; i++) {
                    fx[i] = fx[i-1] + delta[i];
                }
                for (int i = 0; i < delta.length; i++) {
                    fx[i] = Math.exp(fx[i]);
                }
                return fx;
            }

            public double[] incrementsFromParameter(double[] x) {
                double[] increments = new double[x.length];
                increments[0] = Math.log(x[0]);
                for (int i = 1; i < x.length; i++) {
                    increments[i] = Math.log(x[i]/x[i - 1]);
                }
                return increments;
            }
        };

        IncrementTransformType(String transformType) {
            this.transformType = transformType;
        }

        public String getTransformType() {return transformType;}

        private String transformType;
        public abstract double getDerivativeOfInverseTransform(double x);
        public abstract double[] parameterFromIncrements(double[] x);
        public abstract double[] incrementsFromParameter(double[] x);

        public static GradientWrtIncrement.IncrementTransformType factory(String match) {
            for (GradientWrtIncrement.IncrementTransformType type : GradientWrtIncrement.IncrementTransformType.values()) {
                if (match.equalsIgnoreCase(type.getTransformType())) {
                    return type;
                }
            }
            return null;
        }
    }

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

    @Override
    public double[] getGradientLogDensity() {
        // The gradient with respect to the variable-scale
        double[] grad = gradient.getGradientLogDensity();

        // The parameter on the scale of the gradient
        double[] gradScaleParameter = type.parameterFromIncrements(incrementParameter.getParameterValues());

        // The gradient with respect to the increments
        double[] incrementGrad = new double[dim];

        // TODO: is this only right for log-transforms?
        incrementGrad[dim - 1] = grad[dim - 1] * gradScaleParameter[dim - 1];
        for (int i = dim - 2; i > -1; i--) {
            incrementGrad[i] = grad[i] * gradScaleParameter[i] + incrementGrad[i + 1];
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
