package dr.inference.model;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Max Tolkoff
 * @author Marc A. Suchard
 */

public class CompoundGradientProvider implements GradientProvider {

    private int dimension = -1;
    private List<GradientProvider> derivativeList;
    
    public CompoundGradientProvider(List<GradientProvider> derivativeList){
        this.derivativeList = derivativeList;
    }

//    @Override
//    public Likelihood getLikelihood() {
//        return null;
//    }
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
//    public void getGradient(final double[] destination, final int offset) {
//        double[] grad = getGradient();
//        System.arraycopy(grad, 0, destination, offset, grad.length);
//    }

    @Override
    public double[] getGradient() {
        int count = 0;
        ArrayList<double[]> derivativeArrayList = new ArrayList<double[]>();

        for (int i = 0; i < derivativeList.size(); i++) {
            derivativeArrayList.add(derivativeList.get(i).getGradient());
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
