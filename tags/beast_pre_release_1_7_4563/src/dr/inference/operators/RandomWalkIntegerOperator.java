package dr.inference.operators;

import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inferencexml.operators.RandomWalkIntegerOperatorParser;
import dr.math.MathUtils;

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
    double logq;
    public double doOperation() {
        logq = 0.0;
        // a random dimension to perturb
        int index = MathUtils.nextInt(parameter.getSize()); // use getSize(), which = getDimension()

        if (parameter instanceof Parameter) {
            int newValue = calculateNewValue(index);
            ((Parameter) parameter).setParameterValue(index, newValue);
            //System.out.println("newValue: "+newValue);
        } else if (parameter instanceof Variable) { // todo this code is improper if we are going to use Variable<Double>

            int newValue = calculateNewValue(index);
            ((Variable<Integer>) parameter).setValue(index, newValue);

        }

        return logq;
    }

    protected int calculateNewValue(int index) {
        // a random non zero integer around old value within windowSize * 2
        int oldValue;
        int upper;
        int lower;

        if (parameter instanceof Parameter) {
            oldValue = (int) ((Parameter) parameter).getParameterValue(index);
            upper = (int) (double) ((Parameter) parameter).getBounds().getUpperLimit(index);
            lower = (int) (double) ((Parameter) parameter).getBounds().getLowerLimit(index);
        } else if (parameter instanceof Variable) { // todo this code is improper if we are going to use Variable<Double>
            oldValue = ((Variable<Integer>) parameter).getValue(index);
            upper = ((Variable<Integer>) parameter).getBounds().getUpperLimit(index);
            lower = ((Variable<Integer>) parameter).getBounds().getLowerLimit(index);
        } else {
            throw new RuntimeException("The parameter (" + parameter.getId() + ") uses invalid class!");
        }

        if (upper == lower) return upper;

        int maxWindowSize = upper - lower;
        if(windowSize> maxWindowSize){
            windowSize = maxWindowSize;
            System.err.println("The maximum window size should be smaller than the total number of possible integer values.");
        }

        int newValue;
        int roll = MathUtils.nextInt(2 * windowSize); // windowSize="1"; roll = {0, 1}
        if (roll >= windowSize) { // roll = 1
            //roll - window is the positive step size
            int step = 1 + (roll - windowSize);
            newValue = oldValue + step;

            if (newValue > upper){
                newValue = 2 * upper - newValue; //reflect down
            }

        } else {  // roll = 0
            newValue = oldValue - 1 - roll;

            if (newValue < lower){
                newValue = 2 * lower - newValue; //reflect up
            }

        }


        //New and seemingly correct (accoding to the running MCMC with uniform prior)
        //calculation of the hastings ratio --CHW
        int newToOldCount = 0;
        int oldToNewCount = 0;
        if(newValue != oldValue){
            oldToNewCount = oldToNewCount +1;
            newToOldCount = newToOldCount +1;
        }
        int temp = oldValue + windowSize;
        if(temp > upper){
            if((2*upper - temp) <= newValue && newValue != upper){
                oldToNewCount = oldToNewCount+1;
            }
        }

        temp = oldValue - windowSize;
        if(temp < lower){
            if((2*lower - temp) >= newValue && newValue != lower){
                oldToNewCount = oldToNewCount+1;
            }
        }

        temp = newValue + windowSize;
        if( temp > upper){
            if((2*upper - temp) <= oldValue && oldValue != upper){
                newToOldCount = newToOldCount+1;
            }
        }

        temp = newValue - windowSize;
        if( temp < lower){
            if((2*lower - temp) >= oldValue && oldValue != lower){
                newToOldCount = newToOldCount+1;
            }
        }


        logq = Math.log(newToOldCount)- Math.log(oldToNewCount);

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
            maxDelta = ((Parameter) parameter).getParameterValue(0) * 2.0;
        } else if (parameter instanceof Variable) {
            maxDelta = ((Variable<Integer>) parameter).getValue(0) * 2.0;
        }
        long ws = Math.round(OperatorUtils.optimizeWindowSize(windowSize, maxDelta * 2.0, prob, targetProb));

        if (prob < getMinimumGoodAcceptanceLevel()) {
            if(ws <= 1){
                return "";
            }
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
}
