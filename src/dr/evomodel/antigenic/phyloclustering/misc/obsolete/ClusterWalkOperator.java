package dr.evomodel.antigenic.phyloclustering.misc.obsolete;


import dr.inference.model.Bounds;
import dr.inference.model.Parameter;
import dr.inference.operators.*;
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
    public final double doOperation() {
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
            newValue = RandomWalkOperator.reflectValue(newValue, lower, upper);
        } else if (newValue < lower || newValue > upper) {
//            throw new OperatorFailedException("proposed value outside boundaries");
            return Double.NEGATIVE_INFINITY;
        }

        parameter.setParameterValue(index, newValue);

        return 0.0;
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
