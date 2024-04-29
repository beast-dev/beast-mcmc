package dr.evolution.coalescent;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;

/**
 * This is interface for an Interval list that wraps a tree and provides
 * a mapping between nodes in the tree and intervals
 * @author jtmccrone
 */
public interface TreeIntervalList extends IntervalList{
//	void setIntervalStartIndices(int intervalCount);
//	void initializeMaps();

    /**
     *
      * @param nodeNumber the node number
     * @return int[] of interval indices for which the node is a member
     */
    int[] getIntervalsForNode(int nodeNumber);

    /**
     *
     * @param interval the index of interval
     * @return int[] of node numbers that make up the interval
     */
    int[] getNodeNumbersForInterval(int interval);

    /**
     * Set the flag that determines whether or not the interval node mapping should be built.
     * In some cases building the mapping comes with overhead that is best to pass up if the mapping is not needed.
     * @param buildIntervalNodeMapping boolean
     */
    void setBuildIntervalNodeMapping(boolean buildIntervalNodeMapping);

    /**
     * Return the node that ends a coalescent interval. Throws error if the interval is not a coalescent interval
     * @param interval the index of the interval in question
     * @return the node that ends the interval if it is a coalescent interval
     */
    NodeRef getCoalescentNode(int interval);

    /**
     * I think this returns an array with an entry for each internal node sorted
     * by node number. Assuming the array was sorted by height.
     * @param unSortedNodeHeightGradient double[]
     * @return
     */
    double[] sortByNodeNumbers(double[] unSortedNodeHeightGradient);

    /**
     * Get the durations of all intervals that end in a coalescent event.
     * @return array of the coalescent interval durations
     */
    double[] getCoalescentIntervals();

    /**
     * gets the tree
     * @return the tree
     */
    Tree getTree();
}
