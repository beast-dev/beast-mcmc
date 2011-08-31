package dr.inference.operators;

import dr.inference.model.Parameter;
import dr.inferencexml.operators.RandomWalkIntegerOperatorParser;
import dr.math.MathUtils;

import java.util.ArrayList;
import java.util.List;

/**
 * A generic random walk operator for use with a multi-dimensional Integer parameters.
 *
 * @author Michael Defoin Platel
 * @version $Id: RandomWalkIntegerOperator.java$
 */
public class RandomWalkIntegerOperator extends SimpleMCMCOperator {

    public RandomWalkIntegerOperator(Parameter parameter, int windowSize, double weight) {
        this.parameter = parameter;
        this.windowSize = windowSize;
        setWeight(weight);
    }


    public RandomWalkIntegerOperator(Parameter parameter, Parameter updateIndex, int windowSize, double weight) {
        this.parameter = parameter;
        this.windowSize = windowSize;
        setWeight(weight);

        updateMap = new ArrayList<Integer>();
        for (int i = 0; i < updateIndex.getDimension(); i++) {
            if (updateIndex.getParameterValue(i) == 1.0)
                updateMap.add(i);
        }
    }

    /**
     * @return the parameter this operator acts on.
     */
    public Parameter getParameter() {
        return parameter;
    }

    public final int getWindowSize() {
        return windowSize;
    }

    /**
     * change the parameter and return the hastings ratio.
     */
    public double doOperation() {

        // a random dimension to perturb
        int index;
        if (updateMap == null)
            index = MathUtils.nextInt(parameter.getDimension());
        else
            index = updateMap.get(MathUtils.nextInt(updateMap.size()));

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
        return "randomWalkInteger(" + parameter.getParameterName() + ")";
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

    public final String getPerformanceSuggestion() {

        double prob = Utils.getAcceptanceProbability(this);
        double targetProb = getTargetAcceptanceProbability();

        double ws = OperatorUtils.optimizeWindowSize(windowSize, parameter.getParameterValue(0) * 2.0, prob, targetProb);

        if (prob < getMinimumGoodAcceptanceLevel()) {
            return "Try decreasing windowSize to about " + ws;
        } else if (prob > getMaximumGoodAcceptanceLevel()) {
            return "Try increasing windowSize to about " + ws;
        } else return "";
    }

    public String toString() {
        return RandomWalkIntegerOperatorParser.RANDOM_WALK_INT_OP + "(" + parameter.getParameterName() + ", " + windowSize + ", " + getWeight() + ")";
    }

    //PRIVATE STUFF

    protected Parameter parameter = null;
    protected int windowSize = 1;
    protected List<Integer> updateMap = null;
}
