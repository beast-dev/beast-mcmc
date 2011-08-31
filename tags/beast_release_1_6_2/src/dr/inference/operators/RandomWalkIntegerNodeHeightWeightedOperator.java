package dr.inference.operators;

import dr.inference.model.Parameter;
import dr.evomodel.tree.TreeModel;
import dr.evolution.tree.NodeRef;
import dr.math.MathUtils;
import dr.inferencexml.operators.RandomWalkIntegerOperatorParser;
import dr.inferencexml.operators.RandomWalkIntegerNodeHeightWeightedOperatorParser;

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

        return 0.0;
    }

    //MCMCOperator INTERFACE
    public String getOperatorName() {
        return "randomWalkIntegerNodeHeightWeighted(" + parameter.getParameterName() + ")";
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
                "(" + parameter.getParameterName() + ", " + windowSize + ", " + getWeight() + ")";
    }
}
