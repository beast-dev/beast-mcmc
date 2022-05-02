package dr.evomodel.treelikelihood.thorneytreelikelihood;

import dr.evomodel.operators.AbstractTreeOperator;
import dr.evomodel.tree.TreeModel;
import dr.math.MathUtils;

public class ConstrainedTreeOperator extends AbstractTreeOperator {
    public ConstrainedTreeOperator(ConstrainedTreeModel tree, double weight, ConstrainableTreeOperator operator) {
        setWeight(weight);
        constrainedTreeModel = tree;
        this.operator = operator;
        if (tree.getInternalNodeCount() == tree.getSubtreeCount()) {
            throw new IllegalArgumentException(getOperatorName() + " is designed to resolve polytomies; however, the "+
                    "constrained tree is fully resolved. Please remove this operator or provide an unresolved "+
                    "constraints tree.");
        }
        subtreeSizes = new double[tree.getSubtreeCount()];
        for (int i = 0; i < tree.getSubtreeCount(); i++) {
            subtreeSizes[i] = (double) tree.getSubtree(i).getInternalNodeCount()-1; // don't choose subtrees with only 1 internal node. There's no topology to sample.
        }
    }

    @Override
    public String getOperatorName() {
        return "Constrained  " + operator.getOperatorName();
    }

    /**
     * Called by operate(), does the actual operation.
     *
     * @return the hastings ratio
     */
    @Override
    public double doOperation() {

        int subtreeIndex = MathUtils.randomChoicePDF(subtreeSizes);
        TreeModel subtree = constrainedTreeModel.getSubtree(subtreeIndex);
        return operator.doOperation(subtree);
    }

    private final ConstrainedTreeModel constrainedTreeModel;
    private final ConstrainableTreeOperator operator;
    private final double[] subtreeSizes;

}
