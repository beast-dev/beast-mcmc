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

        int[] shuffledInd = MathUtils.shuffled(parameter.getSize()); // use getSize(), which = getDimension()
        int index;
        for(int i = 0; i < sampleSize; i++){
            index = shuffledInd[i];

            int newValue = calculateNewValue(index);
            parameter.setValue(index, newValue);
        }
        return 0.0;
    }

}
