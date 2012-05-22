package dr.inference.operators;

import dr.inference.model.Parameter;
import dr.inferencexml.operators.RandomWalkIntegerNodeHeightWeightedOperatorParser;
import dr.math.MathUtils;

/**
 * @author Chieh-Hsi Wu
 *
 * The probability an internal node is picked to have its state changed depends on the node height.
 */
public class RandomWalkIntegerNodeHeightWeightedOperator extends RandomWalkIntegerOperator{

    private Parameter internalNodeHeights;

    public RandomWalkIntegerNodeHeightWeightedOperator(
            Parameter parameter, int windowSize, double weight, Parameter internalNodeHeights){
        super(parameter, windowSize, weight);
        this.internalNodeHeights = internalNodeHeights;
    }

    public double doOperation() {

        // a random dimension to perturb
        int index = MathUtils.randomChoicePDF(internalNodeHeights.getParameterValues());
      
        int newValue = calculateNewValue(index);
        parameter.setValue(index, newValue);

        return 0.0;
    }

    //MCMCOperator INTERFACE
    public String getOperatorName() {
        return "randomWalkIntegerNodeHeightWeighted(" + parameter.getId() + ")";
    }


    public double getTargetAcceptanceProbability() {
        return 0.234;
    }

    public double getMinimumAcceptanceLevel() {
        return 0.1;
    }

    public double getMaximumAcceptanceLevel() {
        return 0.4;
    }

    public double getMinimumGoodAcceptanceLevel() {
        return 0.20;
    }

    public double getMaximumGoodAcceptanceLevel() {
        return 0.30;
    }
   
    public String toString() {
        return RandomWalkIntegerNodeHeightWeightedOperatorParser.RANDOM_WALK_INT_NODE_HEIGHT_WGT_OP +
                "(" + parameter.getId() + ", " + windowSize + ", " + getWeight() + ")";
    }
}
