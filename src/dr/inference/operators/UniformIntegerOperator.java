package dr.inference.operators;

import dr.inference.model.Parameter;
import dr.inferencexml.operators.UniformIntegerOperatorParser;
import dr.math.MathUtils;

/**
 * A generic uniform sampler/operator for use with a multi-dimensional parameter.
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: UniformOperator.java,v 1.16 2005/06/14 10:40:34 rambaut Exp $
 */
public class UniformIntegerOperator extends SimpleMCMCOperator {

    public UniformIntegerOperator(Parameter parameter, int lower, int upper, double weight) {
        this.parameter = parameter;
        this.lower = lower;
        this.upper = upper;
        setWeight(weight);
    }

    /**
     * @return the parameter this operator acts on.
     */
    public Parameter getParameter() {
        return parameter;
    }

    /**
     * change the parameter and return the hastings ratio.
     */
    public final double doOperation() {

        int index = MathUtils.nextInt(parameter.getDimension());
        int newValue = MathUtils.nextInt(upper - lower) + lower;

        parameter.setParameterValue(index, newValue);

        return 0.0;
    }

    //MCMCOperator INTERFACE
    public final String getOperatorName() {
        return "uniformInteger(" + parameter.getParameterName() + ")";
    }

    public final void optimize(double targetProb) {

        throw new RuntimeException("This operator cannot be optimized!");
    }

    public boolean isOptimizing() {
        return false;
    }

    public void setOptimizing(boolean opt) {
        throw new RuntimeException("This operator cannot be optimized!");
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

    public String getPerformanceSuggestion() {
        if (Utils.getAcceptanceProbability(this) < getMinimumAcceptanceLevel()) {
            return "";
        } else if (Utils.getAcceptanceProbability(this) > getMaximumAcceptanceLevel()) {
            return "";
        } else {
            return "";
        }
    }

    public String toString() {
        return UniformIntegerOperatorParser.UNIFORM_INTEGER_OPERATOR + "(" + parameter.getParameterName() + ")";
    }

    //PRIVATE STUFF

    private Parameter parameter = null;
    private int upper;
    private int lower;
}
