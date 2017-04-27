package dr.inference.model;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Max Tolkoff
 * @author Marc A. Suchard
 */

public class CompoundGradient implements GradientWrtParameterProvider {

    private final int dimension;
    private final List<GradientWrtParameterProvider> derivativeList;
    private final Likelihood likelihood;
    private final Parameter parameter;

    public CompoundGradient(List<GradientWrtParameterProvider> derivativeList) {

        this.derivativeList = derivativeList;

        if (derivativeList.size() == 1) {
            likelihood = derivativeList.get(0).getLikelihood();
            parameter = derivativeList.get(0).getParameter();
            dimension = parameter.getDimension();
        } else {
            List<Likelihood> likelihoodList = new ArrayList<Likelihood>();
            CompoundParameter compoundParameter = new CompoundParameter("hmc");

            int dim = 0;
            for (GradientWrtParameterProvider grad : derivativeList) {
                likelihoodList.add(grad.getLikelihood());
                Parameter p = grad.getParameter();
                compoundParameter.addParameter(p);

                dim += p.getDimension();
            }

            likelihood = new CompoundLikelihood(likelihoodList);
            parameter = compoundParameter;
            dimension = dim;
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

//    @Override
//    public void getGradientLogDensity(final double[] destination, final int offset) {
//        double[] grad = getGradientLogDensity();
//        System.arraycopy(grad, 0, destination, offset, grad.length);
//    }

    @Override
    public double[] getGradientLogDensity() {

        double[] result = new double[dimension];

        int offset = 0;
        for (GradientWrtParameterProvider grad : derivativeList) {
            double[] tmp = grad.getGradientLogDensity();
            System.arraycopy(tmp, 0, result, offset, grad.getDimension());
            offset += grad.getDimension();
        }

        return result;
    }
}
