package dr.evomodel.treelikelihood.thorneytreelikelihood;

import dr.evolution.tree.Tree;
import dr.inference.model.Statistic;

/**
 * This is statistic that returns the subtree root heights from a
 * constrained tree. It can be used to log the tmrca of the clades in the
 * constraints tree.
 */
public class SubtreeRootHeightStatistic extends Statistic.Abstract {
    private ConstrainedTreeModel tree;
    public SubtreeRootHeightStatistic(ConstrainedTreeModel tree) {
        super("subtreetmrca");
        this.tree = tree;
    }
    /**
     * @return the number of dimensions that this statistic has.
     */
    @Override
    public int getDimension() {
        return tree.getSubtreeCount();
    }

    /**
     * @param dim the dimension to return value of
     * @return the statistic's scalar value in the given dimension
     */
    @Override
    public double getStatisticValue(int dim) {
        Tree subtree = tree.getSubtree(dim);
        return subtree.getNodeHeight(subtree.getRoot());
    }
}