package dr.inference.model;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Max Tolkoff
 * @author Marc A. Suchard
 */
public class SumDerivative implements GradientWrtParameterProvider {

    private final int dimension;
    private final Likelihood likelihood;
    private final Parameter parameter;

    List<GradientWrtParameterProvider> derivativeList;

    public SumDerivative(List<GradientWrtParameterProvider> derivativeList){

        // TODO Check that parameters are the same

        this.derivativeList = derivativeList;

        GradientWrtParameterProvider first = derivativeList.get(1);
        dimension = first.getDimension();
        parameter = first.getParameter();

        if (derivativeList.size() == 1) {
            likelihood = first.getLikelihood();
        } else {
            List<Likelihood> likelihoodList = new ArrayList<Likelihood>();

            for (GradientWrtParameterProvider grad : derivativeList) {
                if (grad.getDimension() != dimension) {
                    throw new RuntimeException("Unequal parameter dimensions");
                }
                likelihoodList.add(grad.getLikelihood());
            }
            likelihood = new CompoundLikelihood(likelihoodList);
        }
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
        int size = derivativeList.size();

        final double[] derivative = derivativeList.get(0).getGradientLogDensity();

        for (int i = 1; i < size; i++) {
            final double[] temp = derivativeList.get(i).getGradientLogDensity();
            for (int j = 0; j < temp.length; j++) {
                derivative[j] += temp[j];
            }
        }

        return derivative;
    }
}
