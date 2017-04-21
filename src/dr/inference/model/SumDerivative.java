package dr.inference.model;

import java.util.ArrayList;

/**
 * @author Max Tolkoff
 */
public class SumDerivative implements PotentialDerivativeInterface {
    ArrayList<PotentialDerivativeInterface> derivativeList;

    public SumDerivative(ArrayList<PotentialDerivativeInterface> derivativeList){
        this.derivativeList = derivativeList;
    }


    @Override
    public double[] getDerivative() {
        int size = derivativeList.size();

        double[] derivative = derivativeList.get(0).getDerivative();
        double[] temp;
        for (int i = 1; i < size; i++) {
            temp = derivativeList.get(i).getDerivative();
            for (int j = 0; j < temp.length; j++) {
                derivative[j] += temp[j];
            }
        }

        return derivative;
    }
}
