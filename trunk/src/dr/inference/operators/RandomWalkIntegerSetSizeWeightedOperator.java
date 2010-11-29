package dr.inference.operators;

import dr.evolution.tree.NodeRef;
import dr.evomodel.tree.MicrosatelliteSamplerTreeModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.inferencexml.operators.RandomWalkIntegerNodeHeightWeightedOperatorParser;
import dr.math.MathUtils;

/**
 * @author Chieh-Hsi Wu
 *
 * The probability that an index of the internal states parameter picked
 * is weighted by the number of plausible states of a node given the states of its children.
 */
public class RandomWalkIntegerSetSizeWeightedOperator extends RandomWalkIntegerOperator{

    private MicrosatelliteSamplerTreeModel msatSampleTreeModel;
    private double[] weights;
    private double baseSetSize;
    public RandomWalkIntegerSetSizeWeightedOperator(
        Parameter parameter,
        int windowSize,
        double weight,
        MicrosatelliteSamplerTreeModel msatSampleTreeModel,
        double baseIntervalSize){
        super(parameter, windowSize, weight);
        this.msatSampleTreeModel = msatSampleTreeModel;
        this.baseSetSize = baseIntervalSize;


    }

    private void computeSampleWeights(){
        TreeModel tree = msatSampleTreeModel.getTreeModel();
        int intNodeCount = tree.getInternalNodeCount();
        int extNodeCount = tree.getExternalNodeCount();
        weights = new double[intNodeCount];
        for(int i = 0 ; i < intNodeCount; i++){
            NodeRef node = tree.getNode(i+extNodeCount);
            int lcState = msatSampleTreeModel.getNodeValue(tree.getChild(node, 0));
            int rcState = msatSampleTreeModel.getNodeValue(tree.getChild(node, 1));
            weights[i] = Math.abs(lcState-rcState)+baseSetSize;

        }
    }

    public double doOperation() {
        computeSampleWeights();
        // a random dimension to perturb
        int index = MathUtils.randomChoicePDF(weights);

        int newValue = calculateNewValue(index);
        parameter.setValue(index, newValue);

        return 0.0;
    }

    //MCMCOperator INTERFACE
    public String getOperatorName() {
        return "randomWalkIntegerSetSizeWeighted(" + parameter.getId() + ")";
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
        return RandomWalkIntegerNodeHeightWeightedOperatorParser.RANDOM_WALK_INT_NODE_HEIGHT_WGT_OP + "(" + parameter.getId() + ", " + windowSize + ", " + getWeight() + ")";
    }

}
