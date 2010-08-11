package dr.evomodel.operators;

import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.inference.operators.*;
import dr.inferencexml.operators.RandomWalkOperatorParser;
import dr.math.MathUtils;

/**
 * @author Marc A. Suchard
 * @author John J. Welch
 */
public class FunkyPriorMixerOperator extends AbstractCoercableOperator {

    public FunkyPriorMixerOperator(TreeModel treeModel, double windowSize, RandomWalkOperator.BoundaryCondition bc,
                                   double weight, CoercionMode mode) {

        super(mode);
        this.treeModel = treeModel;
        this.parameter = treeModel.getRootHeightParameter(); // TODO Temporary hack so everything compiles
//        this.parameter = parameter;
        this.windowSize = windowSize;
        this.condition = bc;

        setWeight(weight);
    }

    public final double doOperation() throws OperatorFailedException {

        // a random dimension to perturb
        int index;
        index = MathUtils.nextInt(parameter.getDimension());

        double newValue = parameter.getParameterValue(index) + ((2.0 * MathUtils.nextDouble() - 1.0) * windowSize);

//        treeModel.setNodeHeight(node, 0.0); // TODO How get bounds and set heights properly
//        treeModel.getNodeHeightLower(node);

        double lower = parameter.getBounds().getLowerLimit(index);
        double upper = parameter.getBounds().getUpperLimit(index);

        while (newValue < lower || newValue > upper) {
            if (newValue < lower) {
                if (condition == RandomWalkOperator.BoundaryCondition.reflecting) {
                    newValue = lower + (lower - newValue);
                } else {
                    throw new OperatorFailedException("proposed value outside boundaries");
                }

            }
            if (newValue > upper) {
                if (condition == RandomWalkOperator.BoundaryCondition.reflecting) {
                    newValue = upper - (newValue - upper);
                } else {
                    throw new OperatorFailedException("proposed value outside boundaries");
                }

            }
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
        return RandomWalkOperatorParser.RANDOM_WALK_OPERATOR + "(" + parameter.getParameterName() + ", " + windowSize + ", " + getWeight() + ")";
    }

    private final TreeModel treeModel;
    private final Parameter parameter;
    private double windowSize;
    private final RandomWalkOperator.BoundaryCondition condition;
}
