package dr.inference.model;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Max Tolkoff
 * @author Marc A. Suchard
 */

public class CompoundGradientProvider implements GradientProvider {

    private int dimension = -1;
    private final List<GradientProvider> derivativeList;
    private final Likelihood likelihood;

    public CompoundGradientProvider(List<GradientProvider> derivativeList, List<Likelihood> likelihoodList){
        this.derivativeList = derivativeList;
        if (likelihoodList.size() == 1) {
            likelihood = likelihoodList.get(0);
        } else {
            likelihood = new CompoundLikelihood(likelihoodList);
        }

    }

    @Override
    public Likelihood getLikelihood() {
        return null;
    }
//
//    @Override
//    public Parameter getParameter() {
//        if (compoundParameter == null) {
//            compoundParameter = new CompoundParameter("hmc");
//            for (GradientProvider grad : gradientList) {
//                compoundParameter.addParameter(grad.getParameter());
//            }
//        }
//
//        return compoundParameter;
//    }

    @Override
    public int getDimension() {
        if (dimension == -1) {
            dimension = 0;
            for (GradientProvider grad : derivativeList) {
                dimension += grad.getDimension();
            }
        }
        return dimension;
    }

//    @Override
//    public void getGradientLogDensity(final double[] destination, final int offset) {
//        double[] grad = getGradientLogDensity();
//        System.arraycopy(grad, 0, destination, offset, grad.length);
//    }

    @Override
    public double[] getGradientLogDensity() {
        int count = 0;
        ArrayList<double[]> derivativeArrayList = new ArrayList<double[]>();

        for (int i = 0; i < derivativeList.size(); i++) {
            derivativeArrayList.add(derivativeList.get(i).getGradientLogDensity());
            count += derivativeArrayList.get(i).length;
        }

        double[] answer = new double[count];
        count = 0;
        for (int i = 0; i < derivativeArrayList.size(); i++) {
            System.arraycopy(derivativeArrayList.get(i), 0, answer, count, derivativeArrayList.get(i).length);
            count += derivativeArrayList.get(i).length;
        }
        return answer;
    }
}
