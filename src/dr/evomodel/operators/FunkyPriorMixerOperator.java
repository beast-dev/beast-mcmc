package dr.evomodel.operators;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.operators.FunkyPriorMixerOperatorParser;
import dr.inference.model.Parameter;
import dr.inference.operators.*;
import dr.inferencexml.operators.RandomWalkOperatorParser;
import dr.math.MathUtils;

/**
 * @author Marc A. Suchard
 * @author John J. Welch
 */
public class FunkyPriorMixerOperator extends SimpleMCMCOperator //AbstractCoercableOperator
        implements GibbsOperator {

    public FunkyPriorMixerOperator(TreeModel treeModel, double windowSize, RandomWalkOperator.BoundaryCondition bc,
                                   double weight, CoercionMode mode) {
      
//        super(mode);
        this.treeModel = treeModel;
//        this.parameter = treeModel.getRootHeightParameter();
        this.parameter = null;
//        this.parameter = parameter;
        this.windowSize = windowSize;
        this.condition = bc;

        setWeight(weight);
    }

    public final double doOperation() throws OperatorFailedException {

//        // a random dimension to perturb
//        int index;
//        index = MathUtils.nextInt(parameter.getDimension());
//
//        double newValue = parameter.getParameterValue(index) + ((2.0 * MathUtils.nextDouble() - 1.0) * windowSize);
//
////        treeModel.setNodeHeight(node, 0.0);
////        treeModel.getNodeHeightLower(node);
//
//        double lower = parameter.getBounds().getLowerLimit(index);
//        double upper = parameter.getBounds().getUpperLimit(index);
//
//        while (newValue < lower || newValue > upper) {
//            if (newValue < lower) {
//                if (condition == RandomWalkOperator.BoundaryCondition.reflecting) {
//                    newValue = lower + (lower - newValue);
//                } else {
//                    throw new OperatorFailedException("proposed value outside boundaries");
//                }
//
//            }
//            if (newValue > upper) {
//                if (condition == RandomWalkOperator.BoundaryCondition.reflecting) {
//                    newValue = upper - (newValue - upper);
//                } else {
//                    throw new OperatorFailedException("proposed value outside boundaries");
//                }
//
//            }
//        }
//
//        parameter.setParameterValue(index, newValue);

        double[] minNodeHeights = new double[treeModel.getNodeCount()];
        recursivelyFindNodeMinHeights(treeModel, treeModel.getRoot(), minNodeHeights);

        logForwardDensity = new Double(0.0);
        logBackwardDensity = new Double(0.0);

        recursivelyDrawNodeHeights(treeModel, treeModel.getRoot(), 0.0, 0.0, minNodeHeights); //,
                //logForwardDensity, logBackwardDensity);

//        System.err.println("logFD = " + logForwardDensity);
//        System.err.println("logBD = " + logBackwardDensity);

        return 0.0; //logBackwardDensity - logForwardDensity;    // TODO Think: Is this a Gibbs operator?
    }


    private double recursivelyFindNodeMinHeights(Tree tree, NodeRef node, double[] minNodeHeights) {

        // Post-order traversal
        
        double minHeight;

        if (tree.isExternal(node))
            minHeight = tree.getNodeHeight(node);
        else {
            double minHeightChild0 = recursivelyFindNodeMinHeights(tree, tree.getChild(node, 0), minNodeHeights);
            double minHeightChild1 = recursivelyFindNodeMinHeights(tree, tree.getChild(node, 1), minNodeHeights);
            minHeight = (minHeightChild0 > minHeightChild1) ? minHeightChild0 : minHeightChild1;
        }

        minNodeHeights[node.getNumber()] = minHeight;
        return minHeight;
    }

    private void recursivelyDrawNodeHeights(TreeModel tree, NodeRef node,
                                            double oldParentHeight, double newParentHeight, double[] minNodeHeights) {//,
//                                            Double logForwardDensity, Double logBackwardDensity) {
        // Pre-order traversal

        if (tree.isExternal(node))
            return;

        double oldNodeHeight = tree.getNodeHeight(node);
        double newNodeHeight = oldNodeHeight;

        if (!tree.isRoot(node)) {


            double minHeight = minNodeHeights[node.getNumber()];

            double oldDiff = oldParentHeight - minHeight;
            double newDiff = newParentHeight - minHeight;

            newNodeHeight = MathUtils.nextDouble() * newDiff + minHeight; // Currently uniform
            logForwardDensity -= Math.log(newDiff);
            logBackwardDensity -= Math.log(oldDiff);

//            System.err.println("inner logFD = " + logForwardDensity);

            tree.setNodeHeight(node, newNodeHeight);            
        }

        recursivelyDrawNodeHeights(tree, tree.getChild(node, 0), oldNodeHeight, newNodeHeight, minNodeHeights); //,
//                logForwardDensity, logBackwardDensity);
        recursivelyDrawNodeHeights(tree, tree.getChild(node, 1), oldNodeHeight, newNodeHeight, minNodeHeights); //,
//                logForwardDensity, logBackwardDensity);

    }

    //MCMCOperator INTERFACE

    public final String getOperatorName() {
        return FunkyPriorMixerOperatorParser.FUNKY_OPERATOR;
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

//        double ws = OperatorUtils.optimizeWindowSize(windowSize, parameter.getParameterValue(0) * 2.0, prob, targetProb);
        double ws = 2;

        if (prob < getMinimumGoodAcceptanceLevel()) {
            return "Try decreasing windowSize to about " + ws;
        } else if (prob > getMaximumGoodAcceptanceLevel()) {
            return "Try increasing windowSize to about " + ws;
        } else return "";
    }

    public String toString() {
        return RandomWalkOperatorParser.RANDOM_WALK_OPERATOR + "(" + parameter.getParameterName() + ", " + windowSize + ", " + getWeight() + ")";
    }

    public int getStepCount() {
        return 0;
    }

    private final TreeModel treeModel;
    private final Parameter parameter;
    private double windowSize;
    private final RandomWalkOperator.BoundaryCondition condition;

    private Double logForwardDensity;
    private Double logBackwardDensity;

    /**
     * @return the number of steps the operator performs in one go.
     */
}
