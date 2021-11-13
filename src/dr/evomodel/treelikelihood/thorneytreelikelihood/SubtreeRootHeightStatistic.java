package dr.evomodel.treelikelihood.thorneytreelikelihood;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evomodel.tree.TreeStatistic;
import dr.inference.model.Statistic;

import java.util.List;

/**
 * This is statistic that returns the subtree root heights from a
 * constrained tree. It can be used to log the tmrca of the clades in the
 * constraints tree.
 */
public class SubtreeRootHeightStatistic extends TreeStatistic {
    private ConstrainedTreeModel tree;
    private int dimension;
    private int subtree;
    private final double mostRecentTipTime;
    private final boolean isBackwards;
    public SubtreeRootHeightStatistic(String name, ConstrainedTreeModel tree, TaxonList taxa,boolean isAbsolute) {
        super(name);

        if(taxa==null){
            subtree=-1;
            dimension=tree.getSubtreeCount();
        }else{
            int[] nodes=new int[taxa.getTaxonCount()];
            for (int i = 0; i < nodes.length; i++) {
                nodes[i]=tree.getExternalNode(tree.getTaxonIndex(taxa.getTaxon(i))).getNumber();
            }
            NodeRef tmrca = TreeUtils.getCommonAncestor(tree,nodes);
            subtree=tree.getSubtreeIndex(tmrca);
            dimension=1;
        }

        if (isAbsolute && Taxon.getMostRecentDate() != null) {
            isBackwards = Taxon.getMostRecentDate().isBackwards();
            mostRecentTipTime = Taxon.getMostRecentDate().getAbsoluteTimeValue();
        } else {
            // give node heights or taxa don't have dates
            mostRecentTipTime = Double.NaN;
            isBackwards = false;
        }
        setTree(tree);
    }


        /**
         * @return the number of dimensions that this statistic has.
         */
    @Override
    public int getDimension() {
        return dimension;
    }

    /**
     * @param dim the dimension to return value of
     * @return the statistic's scalar value in the given dimension
     */
    @Override
    public double getStatisticValue(int dim) {
        Tree stree;
        if(subtree==-1){
            stree = tree.getSubtree(dim);
        }else{
            stree = tree.getSubtree(subtree);
        }

        if (!Double.isNaN(mostRecentTipTime)) {
            if (isBackwards) {
                return mostRecentTipTime + stree.getNodeHeight(stree.getRoot());
            } else {
                return mostRecentTipTime - stree.getNodeHeight(stree.getRoot());
            }
        } else {
            return stree.getNodeHeight(stree.getRoot());
        }
    }

    @Override
    public void setTree(Tree tree) {
        this.tree=(ConstrainedTreeModel) tree;
    }

    @Override
    public Tree getTree() {
        return this.tree;
    }
}