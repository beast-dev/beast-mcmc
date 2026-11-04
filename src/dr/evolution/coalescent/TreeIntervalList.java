/*
 * TreeIntervalList.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 *
 */

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

    // This option is set in the constructor as final
//    /**
//     * Set the flag that determines whether or not the interval node mapping should be built.
//     * In some cases building the mapping comes with overhead that is best to pass up if the mapping is not needed.
//     * @param buildIntervalNodeMapping boolean
//     */
//    void setBuildIntervalNodeMapping(boolean buildIntervalNodeMapping);
    boolean isBuildIntervalNodeMapping();

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
