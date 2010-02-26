package dr.inference.operators;

import dr.inference.model.Parameter;
import dr.math.MathUtils;

/**
 * @author Chieh-Hsi Wu
 * Operator that performs independent random walks on k elements in a parameter
 */
public class MultipleRandomWalkIntegerOperator extends RandomWalkIntegerOperator{
    private int sampleSize = 1;
    public MultipleRandomWalkIntegerOperator(Parameter parameter, int windowSize, int sampleSize, double weight) {
        super(parameter,windowSize, weight);
        this.sampleSize = sampleSize;
    }

    public double doOperation() {

        int[] shuffledInd = MathUtils.shuffled(parameter.getDimension());
        int index;
        for(int i = 0; i < sampleSize; i++){

            index = shuffledInd[i];

            // a random non zero integer around old value within windowSize * 2
            int oldValue = (int) parameter.getParameterValue(index);
            int newValue;
            int roll = MathUtils.nextInt(2 * windowSize);
            if (roll >= windowSize) {
                newValue = oldValue + 1 + roll - windowSize;

                if (newValue > parameter.getBounds().getUpperLimit(index))
                    newValue = 2 * (int)(double)parameter.getBounds().getUpperLimit(index) - newValue;
            } else {
                newValue = oldValue - 1 - roll;

                if (newValue < parameter.getBounds().getLowerLimit(index))
                    newValue = 2 * (int)(double)parameter.getBounds().getLowerLimit(index) - newValue;
            }

            parameter.setParameterValue(index, newValue);
        }
        return 0.0;
    }

}
