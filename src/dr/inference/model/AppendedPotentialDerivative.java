package dr.inference.model;

import java.util.ArrayList;

/**
 * @author Max Tolkoff
 */
public class AppendedPotentialDerivative implements PotentialDerivativeInterface {
    ArrayList<PotentialDerivativeInterface> derivativeList;


    public AppendedPotentialDerivative(ArrayList<PotentialDerivativeInterface> derivativeList){
        this.derivativeList = derivativeList;
    }

    @Override
    public double[] getDerivative() {
        int count = 0;
        ArrayList<double[]> derivativeArrayList = new ArrayList<double[]>();

        for (int i = 0; i < derivativeList.size(); i++) {
            derivativeArrayList.add(derivativeList.get(i).getDerivative());
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
