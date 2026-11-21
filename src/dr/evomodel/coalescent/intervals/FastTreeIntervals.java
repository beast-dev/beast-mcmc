/*
 * FastTreeIntervals.java
 *
 * Copyright © 2002-2026, the BEAST Development Team.
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
 */

package dr.evomodel.coalescent.intervals;

import dr.evolution.coalescent.IntervalType;
import dr.evolution.coalescent.TreeIntervalList;
import dr.evolution.tree.MutableTreeModel;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Units;
import dr.evomodel.tree.TreeChangedEvent;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Variable;

import java.util.Arrays;

/**
 * The coalescent intervals of a tree, kept up to date as the tree changes without recalculating them all.
 *
 * Every node of the tree is an event: a tip is a sample and an internal node a coalescence. The events are
 * held in order of time (node height). Interval i lies between event i and event i + 1, so there is one
 * fewer interval than events, interval i starts at the time of event i, and its type and its node are
 * those of the event that ends it, event i + 1. The lineage count of an interval is the number of lineages
 * in the tree during it.
 *
 * The tree's change events say what has to be done before the intervals are next asked for:
 * <ul>
 * <li>the height of a node changed: that node's event is moved to its new place in the order;</li>
 * <li>many, or all, heights changed: every height is read again and the events are put back in order;</li>
 * <li>anything else that changes the nodes: the events are made again from the tree.</li>
 * </ul>
 * A change of topology alone moves no event.
 *
 * Nodes are assumed to be numbered with the tips first, as they are in every tree model.
 *
 * @author JT
 * @author Andrew Rambaut
 */
public class FastTreeIntervals extends AbstractModel implements Units, TreeIntervalList {

    /**
     * What has to be done to the events before they describe the tree, in increasing order of the work
     * involved. A later one takes in the earlier ones.
     */
    private enum Pending {
        NOTHING,
        CHANGED_NODES,
        ALL_TIMES,
        REBUILD
    }

    private final Tree tree;
    private final int intervalCount;

    protected Events events;
    private Events storedEvents;

    private Pending pending;
    private Pending storedPending;

    // the nodes whose heights have changed since the events were last brought up to date
    private int[] changedNodes;
    private int changedNodeCount;
    private int[] storedChangedNodes;
    private int storedChangedNodeCount;

    private Type units = Type.GENERATIONS;

    public FastTreeIntervals(MutableTreeModel tree) {
        this("FastTreeIntervals", tree);
    }

    public FastTreeIntervals(String name, MutableTreeModel tree) {
        // the state is only stored when the tree changes (see AbstractModel)
        super(name, true);

        this.tree = tree;

        final int eventCount = tree.getNodeCount();
        intervalCount = eventCount - 1;

        events = new Events(eventCount, tree.getExternalNodeCount());
        storedEvents = new Events(eventCount, tree.getExternalNodeCount());

        // Beyond this many changed nodes it is cheaper to read all the heights again than to move the
        // events one at a time.
        final int maxChangedNodes = Math.max(8, eventCount / 16);
        changedNodes = new int[maxChangedNodes];
        storedChangedNodes = new int[maxChangedNodes];

        pending = Pending.REBUILD;
        calculateIntervals();

        addModel(tree);
    }

    /**
     * The events will be made again from the tree.
     */
    public void makeDirty() {
        pending = Pending.REBUILD;
        changedNodeCount = 0;
    }

    /**
     * @return whether the events describe the tree as it is now
     */
    protected final boolean isIntervalsKnown() {
        return pending == Pending.NOTHING;
    }

    /**
     * Brings the events up to date with the tree. Every getter does this first, if it is needed.
     */
    @Override
    public void calculateIntervals() {
        switch (pending) {
            case REBUILD:
                events.rebuild(tree);
                break;
            case ALL_TIMES:
                events.updateAllTimes(tree);
                break;
            case CHANGED_NODES:
                for (int i = 0; i < changedNodeCount; i++) {
                    final int node = changedNodes[i];
                    events.updateNode(node, tree.getNodeHeight(tree.getNode(node)));
                }
                break;
            case NOTHING:
                break;
        }
        pending = Pending.NOTHING;
        changedNodeCount = 0;
    }

    private Events currentEvents() {
        if (pending != Pending.NOTHING) {
            calculateIntervals();
        }
        return events;
    }

    // **************************************************************
    // IntervalList IMPLEMENTATION
    // **************************************************************

    @Override
    public int getIntervalCount() {
        currentEvents();
        return intervalCount;
    }

    @Override
    public int getSampleCount() {
        return tree.getTaxonCount();
    }

    /**
     * @return the time of the first event
     */
    @Override
    public double getStartTime() {
        return currentEvents().getTime(0);
    }

    @Override
    public double getInterval(int i) {
        return currentEvents().getInterval(i + 1);
    }

    /**
     * @return the time of the start of interval i, which is the time of event i. The time of the end of
     * the last interval is given by i = getIntervalCount().
     */
    @Override
    public double getIntervalTime(int i) {
        return currentEvents().getTime(i);
    }

    @Override
    public int getLineageCount(int i) {
        return currentEvents().getLineageCount(i + 1);
    }

    @Override
    public int getCoalescentEvents(int i) {
        final Events events = currentEvents();
        if (i < intervalCount - 1) {
            return events.getLineageCount(i + 1) - events.getLineageCount(i + 2);
        } else {
            return events.getLineageCount(i + 1) - 1;
        }
    }

    /**
     * @return the type of the event that ends interval i
     */
    @Override
    public IntervalType getIntervalType(int i) {
        return currentEvents().getType(i + 1);
    }

    /**
     * @return the time of the last event (the height of the root). This is the duration of the intervals
     * only when the first event is at time zero.
     */
    @Override
    public double getTotalDuration() {
        return currentEvents().getTime(intervalCount);
    }

    @Override
    public boolean isBinaryCoalescent() {
        return true;
    }

    @Override
    public boolean isCoalescentOnly() {
        return true;
    }

    // **************************************************************
    // TreeIntervalList IMPLEMENTATION
    // **************************************************************

    /**
     * @param nodeNumber the node number
     * @return the position of the node's event in the order of time. The node ends the interval before
     * this index (if there is one) and starts the interval with this index (if there is one).
     */
    public int getIntervalIndexForNode(int nodeNumber) {
        return currentEvents().getNodePosition(nodeNumber);
    }

    /**
     * The first node starts interval 0, each later node ends one interval and starts the next, and the
     * last node (the root) ends the last interval:
     * <pre>
     * Nodes:            4 ------- 3 ------- 2 ------- 1 ------- 0
     * Intervals:             3         2         1         0
     * IntervalsForNode: [3]     [2,3]     [1,2]     [0,1]      [0]
     * </pre>
     *
     * @param nodeNumber the node number
     * @return the indices of the intervals that the node ends and starts: one index for the first and
     * last nodes, otherwise two.
     */
    @Override
    public int[] getIntervalsForNode(int nodeNumber) {
        final int position = currentEvents().getNodePosition(nodeNumber);
        if (position == 0) {
            return new int[]{0};
        } else if (position == intervalCount) {
            return new int[]{position - 1};
        } else {
            return new int[]{position - 1, position};
        }
    }

    /**
     * @return the numbers of the nodes that start and end interval i
     */
    @Override
    public int[] getNodeNumbersForInterval(int i) {
        final Events events = currentEvents();
        return new int[]{events.getNode(i), events.getNode(i + 1)};
    }

    @Override
    public boolean isBuildIntervalNodeMapping() {
        return true;
    }

    @Override
    public NodeRef getCoalescentNode(int interval) {
        final Events events = currentEvents();
        if (events.getType(interval + 1) != IntervalType.COALESCENT) {
            throw new IllegalArgumentException("interval is not a coalescent interval");
        }
        return tree.getNode(events.getNode(interval + 1));
    }

    public NodeRef getSamplingNode(int interval) {
        final Events events = currentEvents();
        if (events.getType(interval + 1) != IntervalType.SAMPLE) {
            throw new IllegalArgumentException("interval is not a sampling interval");
        }
        return tree.getNode(events.getNode(interval + 1));
    }

    /**
     * @param byIntervalOrder a value for each internal node, in the order of their coalescent events
     * @return the same values in the order of the internal nodes' numbers
     */
    @Override
    public double[] sortByNodeNumbers(double[] byIntervalOrder) {
        final Events events = currentEvents();
        final int tipCount = tree.getExternalNodeCount();
        double[] sortedValues = new double[byIntervalOrder.length];
        int next = 0;
        for (int i = 0; i < events.size(); i++) {
            final int internalNode = events.getNode(i) - tipCount;
            if (internalNode >= 0 && internalNode < sortedValues.length) {
                sortedValues[internalNode] = byIntervalOrder[next];
                next++;
            }
        }
        return sortedValues;
    }

    @Override
    public double[] getCoalescentIntervals() {
        final Events events = currentEvents();
        double[] coalescentIntervals = new double[tree.getInternalNodeCount()];
        int next = 0;
        for (int i = 0; i < intervalCount; i++) {
            if (events.getType(i + 1) == IntervalType.COALESCENT) {
                coalescentIntervals[next] = events.getInterval(i + 1);
                next++;
            }
        }
        return coalescentIntervals;
    }

    @Override
    public Tree getTree() {
        return tree;
    }

    // **************************************************************
    // Units IMPLEMENTATION
    // **************************************************************

    public final Type getUnits() {
        return units;
    }

    public final void setUnits(Type units) {
        this.units = units;
    }

    // **************************************************************
    // Model IMPLEMENTATION
    // **************************************************************

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model != tree) {
            return;
        }

        if (!(tree instanceof TreeModel)) {
            // nothing is known about the events of other kinds of tree
            makeDirty();
        } else if (object instanceof TreeChangedEvent) {
            final TreeChangedEvent event = (TreeChangedEvent) object;
            if (event.isNodeChanged()) {
                // A node event is sent when a node is added to or removed from a branch, or when its
                // height, rate or trait changes. Only a change of height moves an event.
                if (event.isHeightChanged()) {
                    nodeHeightChanged(event.getNode().getNumber());
                }
            } else if (event.isTreeChanged()) {
                if (event.isNodeOrderChanged()) {
                    // the whole tree has changed (this is what an empirical tree distribution sends)
                    makeDirty();
                } else {
                    // The heights have changed but not the nodes or the topology. The order of the
                    // events in time can still have changed (heights scaled past a tip, or set from
                    // height ratios).
                    allHeightsChanged();
                }
            }
        }

        fireModelChanged(object);
    }

    private void nodeHeightChanged(int node) {
        if (pending == Pending.ALL_TIMES || pending == Pending.REBUILD) {
            // this node's height is going to be read again with all the others
            return;
        }
        if (changedNodeCount == changedNodes.length) {
            allHeightsChanged();
        } else {
            changedNodes[changedNodeCount] = node;
            changedNodeCount++;
            pending = Pending.CHANGED_NODES;
        }
    }

    private void allHeightsChanged() {
        if (pending != Pending.REBUILD) {
            pending = Pending.ALL_TIMES;
            changedNodeCount = 0;
        }
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int i, Variable.ChangeType changeType) {
        // there are no variables
    }

    /**
     * Called before the first change to the tree in a step, not at the start of every step.
     */
    @Override
    protected void storeState() {
        storedEvents.copyEvents(events);
        storedPending = pending;
        System.arraycopy(changedNodes, 0, storedChangedNodes, 0, changedNodeCount);
        storedChangedNodeCount = changedNodeCount;
    }

    @Override
    protected void restoreState() {
        Events tmpEvents = storedEvents;
        storedEvents = events;
        events = tmpEvents;

        int[] tmpChangedNodes = storedChangedNodes;
        storedChangedNodes = changedNodes;
        changedNodes = tmpChangedNodes;
        changedNodeCount = storedChangedNodeCount;

        pending = storedPending;
    }

    @Override
    protected void acceptState() {
        // nothing to do
    }

    /**
     * The events of a tree in order of time, with the number of lineages before each.
     */
    protected static class Events {

        private final int eventCount;
        private final int tipCount;

        // the node of each event, in order of time
        private final int[] nodes;
        // the position of each node's event in that order
        private final int[] positions;
        // the time of each event
        private final double[] times;
        // the number of lineages in the interval that ends with each event (0 for the first event)
        private final int[] lineageCounts;

        public Events(int eventCount, int tipCount) {
            this.eventCount = eventCount;
            this.tipCount = tipCount;
            nodes = new int[eventCount];
            positions = new int[eventCount];
            times = new double[eventCount];
            lineageCounts = new int[eventCount];
        }

        public int size() {
            return eventCount;
        }

        /**
         * @return the number of the node of event i
         */
        public int getNode(int i) {
            return nodes[i];
        }

        /**
         * @return the position of the node's event in the order of time
         */
        public int getNodePosition(int node) {
            return positions[node];
        }

        public double getTime(int i) {
            return times[i];
        }

        public IntervalType getType(int i) {
            return isSample(nodes[i]) ? IntervalType.SAMPLE : IntervalType.COALESCENT;
        }

        /**
         * @return the number of lineages in the interval that ends with event i
         */
        public int getLineageCount(int i) {
            return lineageCounts[i];
        }

        /**
         * @return the length of the interval that ends with event i (-1 for the first event, which ends
         * no interval)
         */
        public double getInterval(int i) {
            return i == 0 ? -1 : times[i] - times[i - 1];
        }

        public void copyEvents(Events source) {
            System.arraycopy(source.nodes, 0, nodes, 0, eventCount);
            System.arraycopy(source.positions, 0, positions, 0, eventCount);
            System.arraycopy(source.times, 0, times, 0, eventCount);
            System.arraycopy(source.lineageCounts, 0, lineageCounts, 0, eventCount);
        }

        private boolean isSample(int node) {
            return node < tipCount;
        }

        /**
         * The number of lineages after an event is this much more than before it.
         */
        private int lineageChange(int node) {
            return isSample(node) ? 1 : -1;
        }

        /**
         * Makes the events from the tree. Events at the same time are put in the order of their node
         * numbers, apart from the root, which is last.
         */
        void rebuild(Tree tree) {
            final double[] heights = new double[eventCount];
            Integer[] order = new Integer[eventCount];
            for (int i = 0; i < eventCount; i++) {
                final NodeRef node = tree.getNode(i);
                if (tree.isExternal(node) != isSample(i)) {
                    throw new IllegalArgumentException("The tips of the tree are not numbered before the internal nodes");
                }
                heights[i] = tree.getNodeHeight(node);
                order[i] = i;
            }
            final int root = tree.getRoot().getNumber();
            // a stable sort, so nodes at the same time stay in the order of their numbers
            Arrays.sort(order, (a, b) -> heights[a] == heights[b] ?
                    Boolean.compare(a == root, b == root) :
                    Double.compare(heights[a], heights[b]));

            if (!isSample(order[0])) {
                throw new IllegalArgumentException("The first event is not a sample event");
            }
            for (int i = 0; i < eventCount; i++) {
                nodes[i] = order[i];
                times[i] = heights[order[i]];
            }
            countLineages();
        }

        /**
         * Reads the height of every node again and puts the events back in order of time. The events are
         * sorted by insertion, which takes a single pass when the order has not changed and keeps the
         * existing order of events at the same time.
         */
        void updateAllTimes(Tree tree) {
            for (int i = 0; i < eventCount; i++) {
                times[i] = tree.getNodeHeight(tree.getNode(nodes[i]));
            }

            for (int i = 1; i < eventCount; i++) {
                final double time = times[i];
                if (times[i - 1] > time) {
                    final int node = nodes[i];
                    int j = i;
                    while (j > 0 && times[j - 1] > time) {
                        times[j] = times[j - 1];
                        nodes[j] = nodes[j - 1];
                        j--;
                    }
                    times[j] = time;
                    nodes[j] = node;
                }
            }
            countLineages();
        }

        /**
         * Sets the lineage counts and the positions of the nodes from the order of the events.
         */
        private void countLineages() {
            int lineages = 0;
            for (int i = 0; i < eventCount; i++) {
                lineageCounts[i] = lineages;
                positions[nodes[i]] = i;
                lineages += lineageChange(nodes[i]);
            }
        }

        /**
         * Moves the event of a node whose height has changed to its place in the order. An event that
         * arrives at the time of others is put after them.
         */
        void updateNode(int node, double newTime) {
            final int from = positions[node];
            final double oldTime = times[from];

            int to = from;
            if (newTime > oldTime) {
                to = findFirstLater(newTime, from + 1, eventCount) - 1;
            } else if (newTime < oldTime) {
                to = findFirstLater(newTime, 0, from);
            }

            // The events between the old and new positions each move one place towards the old one.
            // The moved event is now on the other side of them, which changes their lineage counts.
            final int change = lineageChange(node);
            if (to > from) {
                for (int i = from; i < to; i++) {
                    setEvent(i, nodes[i + 1], times[i + 1], lineageCounts[i + 1] - change);
                }
            } else {
                for (int i = from; i > to; i--) {
                    setEvent(i, nodes[i - 1], times[i - 1], lineageCounts[i - 1] + change);
                }
            }
            setEvent(to, node, newTime,
                    to == 0 ? 0 : lineageCounts[to - 1] + lineageChange(nodes[to - 1]));
        }

        private void setEvent(int i, int node, double time, int lineageCount) {
            nodes[i] = node;
            positions[node] = i;
            times[i] = time;
            lineageCounts[i] = lineageCount;
        }

        /**
         * @return the position of the first event from start (inclusive) to end (exclusive) with a time
         * later than the one given, or end if there is none. A binary search.
         */
        private int findFirstLater(double time, int start, int end) {
            int low = start;
            int high = end;
            while (low < high) {
                final int mid = (low + high) >>> 1;
                if (times[mid] <= time) {
                    low = mid + 1;
                } else {
                    high = mid;
                }
            }
            return low;
        }
    }
}
