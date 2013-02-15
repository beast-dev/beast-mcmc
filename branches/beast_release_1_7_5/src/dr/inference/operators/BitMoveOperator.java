package dr.inference.operators;

import dr.inference.model.Parameter;
import dr.math.MathUtils;

import java.util.ArrayList;
import java.util.List;

/**
 * A generic operator that moves k 1 bits to k zero locations.
 *
 * @author Alexei Drummond
 * @version $Id$
 */
public class BitMoveOperator extends SimpleMCMCOperator {

    public BitMoveOperator(Parameter bitsParameter, Parameter valuesParameter, int numBitsToMove, double weight) {
        this.bitsParameter = bitsParameter;
        this.valuesParameter = valuesParameter;

        if (valuesParameter != null && bitsParameter.getDimension() != valuesParameter.getDimension()) {
            throw new IllegalArgumentException("bits parameter must be same length as values parameter");
        }

        this.numBitsToMove = numBitsToMove;
        setWeight(weight);
    }

    /**
     * Pick a random k ones in the vector and move them to a random k zero positions.
     */
    public final double doOperation() throws OperatorFailedException {

        final int dim = bitsParameter.getDimension();
        List<Integer> ones = new ArrayList<Integer>();
        List<Integer> zeros = new ArrayList<Integer>();

        for (int i = 0; i < dim; i++) {
            if (bitsParameter.getParameterValue(i) == 1.0) {
                ones.add(i);
            } else {
                zeros.add(i);
            }
        }

        if (ones.size() >= numBitsToMove && zeros.size() >= numBitsToMove) {

            for (int i = 0; i < numBitsToMove; i++) {

                int myOne = ones.remove(MathUtils.nextInt(ones.size()));
                int myZero = zeros.remove(MathUtils.nextInt(zeros.size()));

                bitsParameter.setParameterValue(myOne, 0.0);
                bitsParameter.setParameterValue(myZero, 1.0);

                if (valuesParameter != null) {
                    double value1 = valuesParameter.getParameterValue(myOne);
                    double value2 = valuesParameter.getParameterValue(myZero);
                    valuesParameter.setParameterValue(myOne, value2);
                    valuesParameter.setParameterValue(myZero, value1);
                }

            }
        } else throw new OperatorFailedException("Not enough bits to move!");

        return 0.0;
    }

    // Interface MCMCOperator
    public final String getOperatorName() {
        StringBuilder builder = new StringBuilder();
        builder.append("bitMove(");
        builder.append(bitsParameter.getParameterName());

        if (valuesParameter != null) {
            builder.append(", ").append(valuesParameter.getParameterName());
        }
        builder.append(", ").append(numBitsToMove).append(")");

        return builder.toString();
    }

    public final String getPerformanceSuggestion() {
        return "no performance suggestion";
    }

    public String toString() {
        return getOperatorName();
    }

    // Private instance variables

    private Parameter bitsParameter = null;
    private Parameter valuesParameter = null;
    private int numBitsToMove = 1;
}
