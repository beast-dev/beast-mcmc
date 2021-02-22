package dr.evomodel.treelikelihood.thorneytreelikelihood;

import dr.evolution.tree.NodeRef;
import dr.evomodel.operators.AbstractAdaptableTreeOperator;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.AdaptationMode;
import dr.math.MathUtils;
import dr.math.Poisson;

import java.util.*;

public class MultiMoveUniformNodeHeightOperator extends AbstractAdaptableTreeOperator {

    private final TreeModel tree;

    /**
     * Constructor
     *
     * @param tree   the tree
     * @param weight the weight
     */
    public MultiMoveUniformNodeHeightOperator(TreeModel tree, double weight, double nonShiftedMean, AdaptationMode mode, double targetAcceptanceProbability){
        super(mode, targetAcceptanceProbability);
        this.tree = tree;
        setWeight(weight);

        this.nonshiftedMean = nonShiftedMean;

    }

    /**
     * Do a subtree leap move.
     *
     * @return the log-transformed hastings ratio
     */
    public double doOperation() {

        List<NodeRef> selected = new ArrayList<>();
        Set<NodeRef> ignore  = new HashSet<>();
        ignore.add(tree.getRoot());
        int operations = getOperationCount();
        int failures = 0;
        while( selected.size()<operations) {

            NodeRef selection =  tree.getInternalNode(MathUtils.nextInt(tree.getInternalNodeCount()));
            if(!ignore.contains(selection)){
                selected.add(selection);
                ignore.add(tree.getParent(selection));
                ignore.add(tree.getChild(selection, 0));
                ignore.add(tree.getChild(selection, 1));
            }else{
                failures++;
            }

            int MAX_FAILURES = 50; //avoiding picking more nodes than reasonable
            if (failures == MAX_FAILURES) {
                return Double.NEGATIVE_INFINITY;
            }
        }
        // Pick a internal nodes avoiding root or children/parent of previously choosen nodes
        double logL = 0.0;
        for (NodeRef node : selected) {
            logL += doMove(node);
        }

        return logL;

    }

    private double doMove(NodeRef node){
        final double upperHeight = tree.getNodeHeight(tree.getParent(node));
        final double lowerHeight = Math.max(
                tree.getNodeHeight(tree.getChild(node, 0)),
                tree.getNodeHeight(tree.getChild(node, 1)));

        final double oldHeight = tree.getNodeHeight(node);

        return doUniform(node, oldHeight, upperHeight, lowerHeight);
    }

    private double doUniform(NodeRef node, double oldValue, double upper, double lower) {
        tree.setNodeHeight(node, (MathUtils.nextDouble() * (upper - lower)) + lower);
        return 0.0;
    }

    public String getOperatorName() {
        return "adaptable uniform(" + tree.getId() + " internal nodes)";
    }

    /**
     * Sets the adaptable parameter value.
     *
     * @param value the value to set the adaptable parameter to
     */
    @Override
    protected void setAdaptableParameterValue(double value) {
        nonshiftedMean = Math.exp(value);
    }

    /**
     * Gets the adaptable parameter value.
     *
     * @returns the value
     */
    @Override
    protected double getAdaptableParameterValue() {
        return Math.log(nonshiftedMean);
    }

    /**
     * @return the underlying tuning parameter value
     */
    @Override
    public double getRawParameter() {
        return nonshiftedMean;
    }
    @Override
    public String getAdaptableParameterName() {
        return "unshiftedMean";
    }

    private int getOperationCount() {
            return Poisson.nextPoisson(this.nonshiftedMean) + 1;
    }

    private double nonshiftedMean;


}
