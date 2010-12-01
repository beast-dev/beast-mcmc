package dr.inference.operators;

import dr.inference.model.Parameter;
import dr.inference.model.Variable;
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

    public RandomWalkIntegerOperator(Variable parameter, int windowSize, double weight) {
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
        return (Parameter) parameter;
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
            index = MathUtils.nextInt(parameter.getSize()); // use getSize(), which = getDimension()
        else
            index = updateMap.get(MathUtils.nextInt(updateMap.size()));

        if (parameter instanceof Parameter) {
            int newValue = calculateNewValue(index);
            ((Parameter) parameter).setParameterValue(index, newValue);
        } else if (parameter instanceof Variable) { // todo this code is improper if we are going to use Variable<Double> 
            if (index > 0) { // if Bounds.Staircase, value of index 0 is fixed to constant 0
                int newValue = calculateNewValue(index);
                ((Variable<Integer>) parameter).setValue(index, newValue);
            }
        }

        return 0.0;
    }

    protected int calculateNewValue(int index) {
        // a random non zero integer around old value within windowSize * 2
        int oldValue;
        int upper;
        int lower;
        if (parameter instanceof Parameter) {
            oldValue = (int) ((Parameter) parameter).getParameterValue(index);
            upper = (int)(double)((Parameter) parameter).getBounds().getUpperLimit(index);
            lower = (int)(double)((Parameter) parameter).getBounds().getLowerLimit(index);
        } else if (parameter instanceof Variable) { // todo this code is improper if we are going to use Variable<Double> 
            oldValue = ((Variable<Integer>) parameter).getValue(index);
            upper = ((Variable<Integer>) parameter).getBounds().getUpperLimit(index);
            lower = ((Variable<Integer>) parameter).getBounds().getLowerLimit(index);
        } else {
            throw new RuntimeException("The parameter (" + parameter.getId() + ") uses invalid class!");
        }
//        System.out.println("index = " + index + ";  oldValue = " + oldValue + ";  upper = " + upper + ";  lower = " + lower);
        if (upper == lower) return upper;

        int newValue;
        int roll = MathUtils.nextInt(2 * windowSize); // windowSize="1"; roll = {0, 1}
        if (roll >= windowSize) { // roll = 1
            newValue = oldValue + 1 + roll - windowSize;

            if (newValue > upper)
                newValue = 2 * upper - newValue;
        } else {  // roll = 0
            newValue = oldValue - 1 - roll;

            if (newValue < lower)
                newValue = 2 * lower - newValue;
        }

        return newValue;
    }

    //MCMCOperator INTERFACE
    public String getOperatorName() {
        return "randomWalkInteger(" + parameter.getId() + ")";
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

        double maxDelta = 0;
        if (parameter instanceof Parameter) {
            maxDelta = ((Parameter) parameter).getParameterValue(0)  * 2.0;
        } else if (parameter instanceof Variable) {
            maxDelta = ((Variable<Integer>) parameter).getValue(0) * 2.0;
        }
        double ws = OperatorUtils.optimizeWindowSize(windowSize, maxDelta * 2.0, prob, targetProb);

        if (prob < getMinimumGoodAcceptanceLevel()) {
            return "Try decreasing windowSize to about " + ws;
        } else if (prob > getMaximumGoodAcceptanceLevel()) {
            return "Try increasing windowSize to about " + ws;
        } else return "";
    }

    public String toString() {
        return RandomWalkIntegerOperatorParser.RANDOM_WALK_INTEGER_OPERATOR + "(" + parameter.getId() + ", " + windowSize + ", " + getWeight() + ")";
    }

    //PRIVATE STUFF

    protected Variable parameter = null;
    protected int windowSize = 1;
    protected List<Integer> updateMap = null;
}
