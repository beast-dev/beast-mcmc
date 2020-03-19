package dr.inference.hmc;

import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.inference.operators.PathDependent;

import java.util.Arrays;

/**
 * @author Marc A. Suchard
 */
public class PathGradient implements GradientWrtParameterProvider, PathDependent {

    private final int dimension;
    private final Likelihood likelihood;
    private final Parameter parameter;

    private final GradientWrtParameterProvider source;
    private final GradientWrtParameterProvider destination;

    private double beta = 1.0;

    public PathGradient(final GradientWrtParameterProvider source,
                        final GradientWrtParameterProvider destination){

        this.source = source;
        this.destination = destination;

        this.dimension = source.getDimension();
        this.parameter = source.getParameter();

//        if (destination.getParameter() != parameter) {
//            throw new RuntimeException("Invalid construction");
//        }
        if (destination.getDimension() != dimension) {
            throw new RuntimeException("Unequal parameter dimensions");
        }
        if (!Arrays.equals(destination.getParameter().getParameterValues(), parameter.getParameterValues())){
            throw new RuntimeException("Unequal parameter values");
        }

        this.likelihood = new Likelihood.Abstract(source.getLikelihood().getModel()) {

            @Override
            protected double calculateLogLikelihood() {

                double likelihood = source.getLikelihood().getLogLikelihood();

                if (beta != 1.0) {
                    likelihood = blend(likelihood, destination.getLikelihood().getLogLikelihood(), beta);
                }

                return likelihood;
            }
        };
    }

    @Override
    public void setPathParameter(double beta) {
        this.beta = beta;
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
        return dimension;
    }

    @Override
    public double[] getGradientLogDensity() {

        final double[] likelihood = source.getGradientLogDensity();

        if (beta != 1.0) {
            
            final double[] destination = this.destination.getGradientLogDensity();

            for (int i = 0; i < likelihood.length; ++i) {
                likelihood[i] = blend(likelihood[i], destination[i], beta);
            }
        }

        return likelihood;
    }

    private static double blend(double source, double destination, double beta) {
        return beta * source + (1.0 - beta) * destination;
    }
}
