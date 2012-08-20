package dr.evomodel.tree;

import dr.evolution.tree.Tree;
import dr.inference.model.Statistic;

/**
 * A statistic that reports the height of a tree
 *
 * @author Alexei Drummond
 * @version $Id: RateStatistic.java,v 1.9 2005/07/11 14:06:25 rambaut Exp $
 */
public class TreeHeightStatistic extends Statistic.Abstract implements TreeStatistic {

    public TreeHeightStatistic(String name, Tree tree) {
        super(name);
        this.tree = tree;
    }

    public void setTree(Tree tree) {
        this.tree = tree;
    }

    public Tree getTree() {
        return tree;
    }

    public int getDimension() {
        return 1;
    }

    /**
     * @return the height of the tree
     */
    public double getStatisticValue(int dim) {

        return tree.getNodeHeight(tree.getRoot());
    }

    private Tree tree = null;
}
