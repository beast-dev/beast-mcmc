package dr.evomodel.antigenic.phyloClustering.misc.obsolete;


import dr.inference.model.Bounds;
import dr.inference.model.Parameter;
import dr.inference.operators.AbstractCoercableOperator;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.OperatorUtils;
import dr.math.MathUtils;

import java.util.ArrayList;
import java.util.List;


public class ClusterWalkOperator extends AbstractCoercableOperator {



    public enum BoundaryCondition {
        reflecting,
        absorbing
    }

    public ClusterWalkOperator(Parameter parameter, double windowSize, BoundaryCondition bc, double weight, CoercionMode mode) {
        this(parameter, null, windowSize, bc, weight, mode);
    }

    public ClusterWalkOperator(Parameter parameter, Parameter updateIndex, double windowSize, BoundaryCondition bc,
                              double weight, CoercionMode mode) {
        this(parameter, updateIndex, windowSize, bc, weight, mode, null, null);
    }

    public ClusterWalkOperator(Parameter parameter, Parameter updateIndex, double windowSize, BoundaryCondition bc,
                              double weight, CoercionMode mode, Double lowerOperatorBound, Double upperOperatorBound) {
        super(mode);
        this.parameter = parameter;
        this.windowSize = windowSize;
        this.condition = bc;

        setWeight(weight);
        if (updateIndex != null) {
            updateMap = new ArrayList<Integer>();
            for (int i = 0; i < updateIndex.getDimension(); i++) {
                if (updateIndex.getParameterValue(i) == 1.0)
                    updateMap.add(i);
            }
        }

        this.lowerOperatorBound = lowerOperatorBound;
        this.upperOperatorBound = upperOperatorBound;
    }

    /**
     * @return the parameter this operator acts on.
     */
    public Parameter getParameter() {
        return parameter;
    }

    public final double getWindowSize() {
        return windowSize;
    }

    /**
     * change the parameter and return the hastings ratio.
     */
    public final double doOperation() throws OperatorFailedException {
    	System.out.println("Walking cluster");

        // a random dimension to perturb
        int index;
        if (updateMap == null) {
            index = MathUtils.nextInt(parameter.getDimension());
        } else {
            index = updateMap.get(MathUtils.nextInt(updateMap.size()));
        }

        // a random point around old value within windowSize * 2
        double draw = (2.0 * MathUtils.nextDouble() - 1.0) * windowSize;
        double newValue = parameter.getParameterValue(index) + draw;

        final Bounds<Double> bounds = parameter.getBounds();
        final double lower = (lowerOperatorBound == null ? bounds.getLowerLimit(index) : Math.max(bounds.getLowerLimit(index), lowerOperatorBound));
        final double upper = (upperOperatorBound == null ? bounds.getUpperLimit(index) : Math.min(bounds.getUpperLimit(index), upperOperatorBound));

        if (condition == BoundaryCondition.reflecting) {
            newValue = reflectValue(newValue, lower, upper);
        } else if (newValue < lower || newValue > upper) {
            throw new OperatorFailedException("proposed value outside boundaries");
        }

        parameter.setParameterValue(index, newValue);

        return 0.0;
    }

    public double reflectValue(double value, double lower, double upper) {

        double newValue = value;

        if (value < lower) {
            if (Double.isInfinite(upper)) {
                // we are only going to reflect once as the upper bound is at infinity...
                newValue = lower + (lower - value);
            } else {
                double remainder = lower - value;

                double widths = Math.floor(remainder / (upper - lower));
                remainder -= (upper - lower) * widths;

                // even reflections
                if (widths % 2 == 0) {
                    newValue = lower + remainder;
                    // odd reflections
                } else {
                    newValue = upper - remainder;
                }
            }
        } else if (value > upper) {
            if (Double.isInfinite(lower)) {
                // we are only going to reflect once as the lower bound is at -infinity...
                newValue = upper - (newValue - upper);
            } else {

                double remainder = value - upper;

                double widths = Math.floor(remainder / (upper - lower));
                remainder -= (upper - lower) * widths;

                // even reflections
                if (widths % 2 == 0) {
                    newValue = upper - remainder;
                    // odd reflections
                } else {
                    newValue = lower + remainder;
                }
            }
        }

        return newValue;
    }

    public double reflectValueLoop(double value, double lower, double upper) {
        double newValue = value;

        while (newValue < lower || newValue > upper) {
            if (newValue < lower) {
                newValue = lower + (lower - newValue);
            }
            if (newValue > upper) {
                newValue = upper - (newValue - upper);

            }
        }

        return newValue;
    }

    //MCMCOperator INTERFACE
    public final String getOperatorName() {
        return parameter.getParameterName();
    }

    public double getCoercableParameter() {
        return Math.log(windowSize);
    }

    public void setCoercableParameter(double value) {
        windowSize = Math.exp(value);
    }

    public double getRawParameter() {
        return windowSize;
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

        double prob = MCMCOperator.Utils.getAcceptanceProbability(this);
        double targetProb = getTargetAcceptanceProbability();

        double ws = OperatorUtils.optimizeWindowSize(windowSize, parameter.getParameterValue(0) * 2.0, prob, targetProb);

        if (prob < getMinimumGoodAcceptanceLevel()) {
            return "Try decreasing windowSize to about " + ws;
        } else if (prob > getMaximumGoodAcceptanceLevel()) {
            return "Try increasing windowSize to about " + ws;
        } else return "";
    }

    public String toString() {
        return ClusterWalkOperatorParser.CLUSTER_WALK_OPERATOR + "(" + parameter.getParameterName() + ", " + windowSize + ", " + getWeight() + ")";
    }

    //PRIVATE STUFF

    private Parameter parameter = null;
    private double windowSize = 0.01;
    private List<Integer> updateMap = null;
    private final BoundaryCondition condition;

    private final Double lowerOperatorBound;
    private final Double upperOperatorBound;
}
