package dr.inference.model;

import java.util.List;

/**
 * @author Max Tolkoff
 */
public class SumDerivative implements GradientProvider {

    private int dimension = -1;

    List<GradientProvider> derivativeList;

    public SumDerivative(List<GradientProvider> derivativeList){
        this.derivativeList = derivativeList;

        for (GradientProvider grad : derivativeList) {
            if (dimension == -1) {
                dimension = grad.getDimension();
            } else {
                if (dimension != grad.getDimension()) {
                    throw new IllegalArgumentException("All gradients must have the same dimension");
                }
            }
        }
    }

    @Override
    public int getDimension() {
        return dimension;
    }

    @Override
    public double[] getGradient() {
        int size = derivativeList.size();

        double[] derivative = derivativeList.get(0).getGradient();
        double[] temp;
        for (int i = 1; i < size; i++) {
            temp = derivativeList.get(i).getGradient();
            for (int j = 0; j < temp.length; j++) {
                derivative[j] += temp[j];
            }
        }

        return derivative;
    }
}
