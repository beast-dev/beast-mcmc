package dr.inference.hmc;

import dr.evomodel.speciation.SpeciationLikelihoodGradient;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
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

            public double[] parameterFromIncrements(double[] x) {
                double[] fx = new double[x.length];
                fx[0] = x[0];
                for (int i = 1; i < x.length; i++) {
                    fx[i] = fx[i-1] + x[i];
                }
                for (int i = 0; i < x.length; i++) {
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

        // Accumulate over vector
        incrementGrad[0] = grad[0] * type.getDerivativeOfInverseTransform(gradScaleParameter[0]);
        for (int i = 1; i < dim; i++) {
            incrementGrad[i] = incrementGrad[i - 1] + grad[i] * type.getDerivativeOfInverseTransform(gradScaleParameter[i]);
        }

        return incrementGrad;
    }

    @Override
    public String getReport() {
        StringBuilder sb = new StringBuilder();

        sb.append("Gradient WRT increments: ").append(new dr.math.matrixAlgebra.Vector(getGradientLogDensity())).append("\n");
        sb.append("Gradient WRT parameters: ").append(new dr.math.matrixAlgebra.Vector(gradient.getGradientLogDensity())).append("\n");
        sb.append("Increments: ").append(new dr.math.matrixAlgebra.Vector(incrementParameter.getParameterValues())).append("\n");
        sb.append("Parameters: ").append(new dr.math.matrixAlgebra.Vector(type.parameterFromIncrements(gradient.getParameter().getParameterValues()))).append("\n");

        return sb.toString();
    }
}
