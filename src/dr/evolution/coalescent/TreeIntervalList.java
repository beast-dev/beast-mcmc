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
 * @author arambaut
 */
public interface TreeIntervalList extends IntervalList {

    /**
     * Has the interval changed since the last tree update?
     * @param i the interval
     */
    boolean hasIntervalChanged(int i);

        /**
         * Returns the type of event at the start of interval.
         */
    IntervalType getIntervalStartType(int i);

    /**
     * Returns the type of event at the end of the interval.
     */
    IntervalType getIntervalEndType(int i);

    /**
     * Returns the time of the start of an interval
     */
    double getIntervalStartTime(int i);

    /**
     * Returns the time of the end of an interval
     */
    double getIntervalEndTime(int i);
    /**
     * Returns the node number at the end of the interval.
     */
    int getIntervalNodeNumber(int i);

    /**
     * Returns the node number at the end of the interval.
     */
    int getIntervalEndNodeNumber(int i);

    /**
     * Returns the node number at the start of the interval.
     */
    int getIntervalStartNodeNumber(int i);

    /**
     * Gets the node ref for the node at the end of the interval
     * @param interval the interval index
     * @return the node ref
     */
    NodeRef getIntervalNode(int interval);

    /**
     * Gets the node ref for the node at the start of the interval
     * @param interval the interval index
     * @return the node ref
     */
    NodeRef getIntervalStartNode(int interval);

    /**
     * Gets the node ref for the node at the end of the interval
     * @param interval the interval index
     * @return the node ref
     */
    NodeRef getIntervalEndNode(int interval);

    /**
     * gets the tree
     * @return the tree
     */
    Tree getTree();
}
