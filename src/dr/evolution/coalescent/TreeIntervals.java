/*
 * TreeIntervals.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
 */

package dr.evolution.coalescent;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.util.HeapSort;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Extracts the intervals from a tree.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 */
public class TreeIntervals implements IntervalList {

    /**
     * Parameterless constructor.
     */
    public TreeIntervals() {
    }


    public TreeIntervals(Tree tree) {
        setTree(tree);
    }


    /**
     * Sets the tree for which intervals are obtained
     */
    public void setTree(Tree tree) {
        this.tree = tree;
        intervalsKnown = false;
    }

    /**
     * Specifies that the intervals are unknown (i.e., the tree has changed).
     */
    public void setIntervalsUnknown() {
        intervalsKnown = false;
    }

    /**
     * Sets the limit for which adjacent events are merged.
     *
     * @param multifurcationLimit A value of 0 means merge addition of leafs (terminal nodes) when possible but
     *                            return each coalescense as a separate event.
     */
    public void setMultifurcationLimit(double multifurcationLimit) {
        this.multifurcationLimit = multifurcationLimit;
        intervalsKnown = false;
    }

    public int getSampleCount() {
        return tree.getExternalNodeCount();
    }

    /**
     * get number of intervals
     */
    public int getIntervalCount() {
        if (!intervalsKnown) {
            calculateIntervals();
        }
        return intervalCount;
    }

    /**
     * Gets an interval.
     */
    public double getInterval(int i) {
        if (!intervalsKnown) {
            calculateIntervals();
        }
        if (i >= intervalCount) throw new IllegalArgumentException();
        return intervals[i];
    }

    /**
     * Returns the number of uncoalesced lineages within this interval.
     * Required for s-coalescents, where new lineages are added as
     * earlier samples are come across.
     */
    public int getLineageCount(int i) {
        if (!intervalsKnown) {
            calculateIntervals();
        }
        if (i >= intervalCount) throw new IllegalArgumentException();
        return lineageCounts[i];
    }

    /**
     * @param interval
     * @return a list of the noderefs representing the lineages in the ith interval.
     */
    public final List getLineages(int interval) {

        if (lineages[interval] == null) {

            List<Object> lines = new ArrayList<Object>();
            for (int i = 0; i <= interval; i++) {
                if (lineagesAdded[i] != null) lines.addAll(lineagesAdded[i]);
                if (lineagesRemoved[i] != null) lines.removeAll(lineagesRemoved[i]);
            }
            lineages[interval] = Collections.unmodifiableList(lines);

        }
        return lineages[interval];
    }

    /**
     * Returns the number coalescent events in an interval
     */
    public int getCoalescentEvents(int i) {
        if (!intervalsKnown) {
            calculateIntervals();
        }
        if (i >= intervalCount) throw new IllegalArgumentException();
        if (i < intervalCount - 1) {
            return lineageCounts[i] - lineageCounts[i + 1];
        } else {
            return lineageCounts[i] - 1;
        }
    }

    /**
     * Returns the type of interval observed.
     */
    public IntervalType getIntervalType(int i) {
        if (!intervalsKnown) {
            calculateIntervals();
        }
        if (i >= intervalCount) throw new IllegalArgumentException();
        int numEvents = getCoalescentEvents(i);

        if (numEvents > 0) return IntervalType.COALESCENT;
        else if (numEvents < 0) return IntervalType.SAMPLE;
        else return IntervalType.NOTHING;
    }

    public NodeRef getCoalescentNode(int interval) {
        if (getIntervalType(interval) == IntervalType.COALESCENT) {
            if (lineagesRemoved[interval] != null) {
                if (lineagesRemoved[interval].size() == 1) {
                    return lineagesRemoved[interval].get(0);
                } else throw new IllegalArgumentException("multiple lineages lost over this interval!");
            } else throw new IllegalArgumentException("Inconsistent: no intervals lost over this interval!");
        } else throw new IllegalArgumentException("Interval " + interval + " is not a coalescent interval.");
    }

    /**
     * get the total height of the genealogy represented by these
     * intervals.
     */
    public double getTotalDuration() {

        if (!intervalsKnown) {
            calculateIntervals();
        }
        double height = 0.0;
        for (int j = 0; j < intervalCount; j++) {
            height += intervals[j];
        }
        return height;
    }

    /**
     * Checks whether this set of coalescent intervals is fully resolved
     * (i.e. whether is has exactly one coalescent event in each
     * subsequent interval)
     */
    public boolean isBinaryCoalescent() {
        if (!intervalsKnown) {
            calculateIntervals();
        }
        for (int i = 0; i < intervalCount; i++) {
            if (getCoalescentEvents(i) > 0) {
                if (getCoalescentEvents(i) != 1) return false;
            }
        }

        return true;
    }

    /**
     * Checks whether this set of coalescent intervals coalescent only
     * (i.e. whether is has exactly one or more coalescent event in each
     * subsequent interval)
     */
    public boolean isCoalescentOnly() {
        if (!intervalsKnown) {
            calculateIntervals();
        }
        for (int i = 0; i < intervalCount; i++) {
            if (getCoalescentEvents(i) < 1) return false;
        }

        return true;
    }

    /**
     * Recalculates all the intervals for the given tree.
     */
    private void calculateIntervals() {

        int nodeCount = tree.getNodeCount();

        times = new double[nodeCount];
        int[] childCounts = new int[nodeCount];

        collectTimes(tree, times, childCounts);

        indices = new int[nodeCount];

        HeapSort.sort(times, indices);

        if (intervals == null || intervals.length != nodeCount) {
            intervals = new double[nodeCount];
            lineageCounts = new int[nodeCount];
            lineagesAdded = new List[nodeCount];
            lineagesRemoved = new List[nodeCount];
            lineages = new List[nodeCount];
        }

        // start is the time of the first tip
        double start = times[indices[0]];
        int numLines = 0;
        int nodeNo = 0;
        intervalCount = 0;
        while (nodeNo < nodeCount) {

            int lineagesRemoved = 0;
            int lineagesAdded = 0;

            double finish = times[indices[nodeNo]];
            double next;

            do {
                final int childIndex = indices[nodeNo];
                final int childCount = childCounts[childIndex];
                // dont use nodeNo from here on in do loop
                nodeNo += 1;
                if (childCount == 0) {
                    addLineage(intervalCount, tree.getNode(childIndex));
                    lineagesAdded += 1;
                } else {
                    lineagesRemoved += (childCount - 1);

                    // record removed lineages
                    final NodeRef parent = tree.getNode(childIndex);
                    //assert childCounts[indices[nodeNo]] == tree.getChildCount(parent);
                    //for (int j = 0; j < lineagesRemoved + 1; j++) {
                    for (int j = 0; j < childCount; j++) {
                        NodeRef child = tree.getChild(parent, j);
                        removeLineage(intervalCount, child);
                    }

                    // record added lineages
                    addLineage(intervalCount, parent);
                    // no mix of removed lineages when 0 th
                    if (multifurcationLimit == 0.0) {
                        break;
                    }
                }

                if (nodeNo < nodeCount) {
                    next = times[indices[nodeNo]];
                } else break;
            } while (Math.abs(next - finish) <= multifurcationLimit);

            if (lineagesAdded > 0) {

                if (intervalCount > 0 || ((finish - start) > multifurcationLimit)) {
                    intervals[intervalCount] = finish - start;
                    lineageCounts[intervalCount] = numLines;
                    intervalCount += 1;
                }

                start = finish;
            }

            // add sample event
            numLines += lineagesAdded;

            if (lineagesRemoved > 0) {

                intervals[intervalCount] = finish - start;
                lineageCounts[intervalCount] = numLines;
                intervalCount += 1;
                start = finish;
            }
            // coalescent event
            numLines -= lineagesRemoved;
        }

        intervalsKnown = true;
    }

    /**
     * Returns the time of the start of an interval
     *
     * @param i which interval
     * @return start time
     */
    public double getIntervalTime(int i) {
        if (!intervalsKnown) {
            calculateIntervals();
        }
        return times[indices[i]];
    }

    private void addLineage(int interval, NodeRef node) {
        if (lineagesAdded[interval] == null) lineagesAdded[interval] = new ArrayList<NodeRef>();
        lineagesAdded[interval].add(node);
    }

    private void removeLineage(int interval, NodeRef node) {
        if (lineagesRemoved[interval] == null) lineagesRemoved[interval] = new ArrayList<NodeRef>();
        lineagesRemoved[interval].add(node);
    }

    /**
     * @return the delta parameter of Pybus et al (Node spread statistic)
     */
    public double getDelta() {

        return IntervalList.Utils.getDelta(this);
    }

    /**
     * extract coalescent times and tip information into array times from tree.
     */
    private static void collectTimes(Tree tree, double[] times, int[] childCounts) {

        for (int i = 0; i < tree.getNodeCount(); i++) {
            NodeRef node = tree.getNode(i);
            times[i] = tree.getNodeHeight(node);
            childCounts[i] = tree.getChildCount(node);
        }
    }

    /**
     * Return the units that this tree is expressed in.
     */
    public final Type getUnits() {
        return tree.getUnits();
    }

    /**
     * Sets the units that this tree is expressed in.
     */
    public final void setUnits(Type units) {
        throw new IllegalArgumentException("Can't set interval's units");
    }

    /**
     * Extra functionality to store and restore values for caching
     */
    public void storeState() {
        if (intervalsKnown) {

            if (storedIntervals == null) {
                storedIntervals = new double[intervals.length];
            }
            if (storedLineageCounts == null) {
                storedLineageCounts = new int[lineageCounts.length];
            }
            if (superStore) {
            if (storedIndices == null) {
                storedIndices = new int[indices.length];
            }
            if (storedTimes == null) {
                storedTimes = new double[times.length];
            }
            }
            System.arraycopy(intervals, 0, storedIntervals, 0, intervals.length);
            System.arraycopy(lineageCounts, 0, storedLineageCounts, 0, lineageCounts.length);
            if (superStore) {
            System.arraycopy(indices, 0, storedIndices, 0, indices.length);
            System.arraycopy(times, 0, storedTimes, 0, times.length);
            }
        }

        storedIntervalsKnown = intervalsKnown;
    }

    public void restoreState() {
        intervalsKnown = storedIntervalsKnown;

        if (intervalsKnown) {
            double[] tmp1 = storedIntervals;
            storedIntervals = intervals;
            intervals = tmp1;

            int[] tmp2 = storedLineageCounts;
            storedLineageCounts = lineageCounts;
            lineageCounts = tmp2;

            if (superStore) {
            tmp1 = storedTimes;
            storedTimes = times;
            times = tmp1;

            tmp2 = storedIndices;
            storedIndices = indices;
            indices = tmp2;
            }
        }
    }

    public String toString() {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < getIntervalCount(); i++) {
            sb.append("[ ");
            sb.append(getInterval(i));
            sb.append(": ");
            sb.append(getIntervalTime(i));
            sb.append(": ");
            sb.append(getLineageCount(i));
            sb.append(" ]");
        }
        return sb.toString();
    }

    private int[] indices;
    private int[] storedIndices;

    private double[] times;
    private double[] storedTimes;

    /**
     * The tree.
     */
    private Tree tree = null;

    /**
     * The widths of the intervals.
     */
    private double[] intervals;
    private double[] storedIntervals;

    /**
     * The number of uncoalesced lineages within a particular interval.
     */
    private int[] lineageCounts;
    private int[] storedLineageCounts;

    /**
     * The lineages in each interval (stored by node ref).
     */
    private List<NodeRef>[] lineagesAdded;
    private List<NodeRef>[] lineagesRemoved;
    private List[] lineages;

    private int intervalCount = 0;

    /**
     * are the intervals known?
     */
    private boolean intervalsKnown = false;
    private boolean storedIntervalsKnown;
	
	private double multifurcationLimit = -1.0;

    private static final boolean superStore = true;
}