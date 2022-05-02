package dr.evomodel.bigfasttree;

import dr.evolution.coalescent.IntervalList;
import dr.evolution.coalescent.IntervalType;
import dr.evolution.tree.NodeRef;
import dr.evolution.util.Units;
import dr.evomodel.tree.TreeChangedEvent;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Variable;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class BigFastTreeIntervals extends AbstractModel implements Units, IntervalList {
    public BigFastTreeIntervals(TreeModel tree) {
        super("bigFastIntervals");
        int maxEventCount = tree.getNodeCount();

        this.tree = tree;
        intervalsKnown = false;

        updatedNodes = new ArrayList<>(maxEventCount);
        storedUpdatedNodes = new ArrayList<>(maxEventCount);

        events = new Events(maxEventCount);
        storedEvents = new Events(maxEventCount);

        dirty = true;
        onlyUpdateTimes = false;
        calculateIntervals();
        addModel(tree);
    }


    public void makeDirty() {
        dirty = true;
    }

    @Override
    public int getIntervalCount() {
        if (!intervalsKnown) {
            calculateIntervals();
        }
        return intervalCount;
    }

    @Override
    public int getSampleCount() {
        return tree.getTaxonCount();
    }


    @Override
    public double getStartTime() {
        if (!intervalsKnown) {
            calculateIntervals();
        }
        return events.getTime(0);
    }

    @Override
    public double getInterval(int i) {
        if (!intervalsKnown) {
            calculateIntervals();
        }
        return events.getInterval(i + 1);
    }

    @Override
    public double getIntervalTime(int i) {
        if (!intervalsKnown) {
            calculateIntervals();
        }
        return events.getTime(i);
    }

    @Override
    public int getLineageCount(int i) {
        if (!intervalsKnown) {
            calculateIntervals();
        }
        return events.getLineageCount(i + 1);
    }

    @Override
    public int getCoalescentEvents(int i) {
        if (!intervalsKnown) {
            calculateIntervals();
        }
        if (i < intervalCount - 1) {
            return events.getLineageCount(i + 1) - events.getLineageCount(i + 2);
        } else {
            return events.getLineageCount(i + 1) - 1;
        }
    }

    @Override
    public IntervalType getIntervalType(int i) {
        if (!intervalsKnown) {
            calculateIntervals();
        }
        return events.getType(i + 1);
    }

    @Override
    public double getTotalDuration() {
        if (!intervalsKnown) {
            calculateIntervals();
        }
        return events.getTime(events.size() - 1);
    }

    @Override
    public boolean isBinaryCoalescent() {
        return true;
    }

    @Override
    public boolean isCoalescentOnly() {
        return true;
    }

    @Override
    public void calculateIntervals() {
        //If dirty we rebuild the evens and sort them using parallel sort
        if (dirty) {
            // Resort nodes by heights
            // will this update the tree nodes?

            NodeRef[] nodes = new NodeRef[tree.getNodeCount()];
            System.arraycopy(tree.getNodes(), 0, nodes, 0, tree.getNodeCount());
            Arrays.parallelSort(nodes, (a, b) -> Double.compare(tree.getNodeHeight(a), tree.getNodeHeight(b)));

            intervalCount = nodes.length - 1;

            double lastTime = tree.getNodeHeight(nodes[0]);
            if (!tree.isExternal(nodes[0])) {
                throw new IllegalArgumentException("The first event is not a sample event");
            }
            events.setEvent(lastTime, IntervalType.SAMPLE, nodes[0].getNumber(), 0,-1, 0);

            int lineages = 1;
            for (int i = 1; i < nodes.length; i++) {
                NodeRef node = nodes[i];
                double time = tree.getNodeHeight(node);
                double interval = time - lastTime;
                IntervalType type = tree.isExternal(node) ? IntervalType.SAMPLE : IntervalType.COALESCENT;
                int lineageCount = lineages;
                events.setEvent(time, type, node.getNumber(), interval, lineageCount, i);
                events.setNodeOrder(node.getNumber(), i);

                if (type == IntervalType.SAMPLE) {
                    lineages++;
                } else {
                    lineages--;
                }
                lastTime = time;
            }
            intervalsKnown = true;

        } else if (onlyUpdateTimes) {
            for (int i = 0; i < events.size(); i++) {
                double newTime = tree.getNodeHeight(tree.getNode(events.getNode(i)));
                events.updateEventTime(newTime, i);
            }
            onlyUpdateTimes = false;
        } else {
            for (int node : updatedNodes) {
                events.updateForChangedNode(node, tree.getNodeHeight(tree.getNode(node)));
            }
        }

        intervalsKnown = true;
        dirty = false;
        updatedNodes = new ArrayList<>();
    }

    private Type units = Type.GENERATIONS;

    public final Type getUnits() {
        return units;
    }

    public final void setUnits(Type units) {
        this.units = units;
    }


    // **************************************************************
    // ModelListener IMPLEMENTATION
    // **************************************************************

    protected void handleModelChangedEvent(Model model, Object object, int index) {

        if (model == tree) {
            if (object instanceof TreeChangedEvent) {
                TreeChangedEvent treeChangedEvent = (TreeChangedEvent) object;
                if (treeChangedEvent.isNodeChanged()) {

                    //  A node event is fired if a node is added to a branch, removed from a branch or its height or
                    // rate changes.
                    if (treeChangedEvent.isHeightChanged()) {
                        NodeRef node = ((TreeChangedEvent) object).getNode();
                        updatedNodes.add(node.getNumber());
                        intervalsKnown = false;
                    }

                } else if (treeChangedEvent.isTreeChanged()) {
                    if (!treeChangedEvent.isNodeOrderChanged()) {
                        onlyUpdateTimes = true;
                    } else {
                        // Full tree events result in a complete updating of the tree likelihood
                        // This event type is now used for EmpiricalTreeDistributions.
//                    System.err.println("Full tree update event - these events currently aren't used\n" +
//                            "so either this is in error or a new feature is using them so remove this message.");
                        makeDirty();
                    }

                }  // Other event types are ignored (probably trait changes).
                //System.err.println("Another tree event has occured (possibly a trait change).");
            }

            fireModelChanged();
        }
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int i, Variable.ChangeType changeType) {

    }

    @Override
    protected void storeState() {
        storedUpdatedNodes = new ArrayList<>();
        storedUpdatedNodes.addAll(updatedNodes);
//        Collections.copy(storedUpdatedNodes, updatedNodes);
        storedIntervalsKnown = intervalsKnown;
        storedEvents.copyEvents(events);
        storedOnlyUpdateTimes = onlyUpdateTimes;
    }

    @Override
    protected void restoreState() {

        List<Integer> tmp2 = storedUpdatedNodes;
        storedUpdatedNodes = updatedNodes;
        updatedNodes = tmp2;
        events.copyEvents(storedEvents);
        intervalsKnown = storedIntervalsKnown;

        onlyUpdateTimes = storedOnlyUpdateTimes;

    }

    @Override
    protected void acceptState() {

    }

    /**
     * A private data class for sorting events on a full recalculation
     * keeps the time, node number and type of event all in one place.
     */
    private class Event implements Comparable<Event> {

        public int compareTo(Event e) {
            double t = e.time;
            if (t < time) {
                return 1;
            } else if (t > time) {
                return -1;
            } else {
                // events are at exact same time so sort by type
                return type.compareTo(e.type);
            }
        }

        Event(double time, IntervalType type, int node, double interval, int lineageCount) {
            this.time = time;
            this.type = type;
            this.node = node;
            this.interval = interval;
            this.lineageCount = lineageCount;
//            this.info = info;
        }

        /* The type of event
         */
        final IntervalType type;

        /**
         * The time of the event
         */
        final double time;
        final double interval;
        final int lineageCount;

        /**
         * Some extra information for the event (e.g., destination of a migration)
         */
//        final int info;

        final int node;
    }

    /**
     * A private classs that wraps an array of events and provides some a higher level api for updated an event when
     * the associated node's height changes.
     */
    private class Events {

        public Events(int numberOfEvents) {
            nodes = new int[numberOfEvents];
            nodeOrder = new int[numberOfEvents];
            lineageCounts = new int[numberOfEvents];
            intervals = new double[numberOfEvents];
            times = new double[numberOfEvents];
            intervalTypes = new IntervalType[numberOfEvents];
            this.numberOfEvents = numberOfEvents;
        }


        public int getNodePosition(int nodeNum) {
            return nodeOrder[nodeNum];
        }

        public IntervalType getType(int i) {
            return intervalTypes[i];
        }

        public double getTime(int i) {
            return times[i];
        }

        public int size() {
            return numberOfEvents;
        }

        public void setEvent(double time, IntervalType type, int node, double interval, int lineageCount, int i) {
            times[i] = time;
            intervalTypes[i] = type;
            nodes[i] = node;
            intervals[i] = interval;
            lineageCounts[i] = lineageCount;
        }

        public void updateForChangedNode(int node, double newTime) {
            int oldPostion = nodeOrder[node];
            double oldHieght = times[oldPostion];

            int newPosition = oldPostion;
            if (newTime > oldHieght) {
                int firstHeigher = findFirstGreater(times, newTime, oldPostion, numberOfEvents - 1);
                newPosition = firstHeigher - 1;

            } else if (newTime < oldHieght) {
                newPosition = findFirstGreater(times, newTime, 0, oldPostion - 1);
            }


            if (newPosition != oldPostion) {
                // Borrowed from SortableArray
                int lineageCount = lineageCounts[oldPostion];
                double interval = intervals[oldPostion];
                double time = times[oldPostion];
                IntervalType type = intervalTypes[oldPostion];

                if (newPosition > oldPostion) {
                    int len = newPosition - oldPostion;
                    // move everything from the src to the dest down one spot (opening up a gap at dest)
                    System.arraycopy(intervals, oldPostion + 1, intervals, oldPostion, len);
                    System.arraycopy(times, oldPostion + 1, times, oldPostion, len);
                    System.arraycopy(intervalTypes, oldPostion + 1, intervalTypes, oldPostion, len);
                    System.arraycopy(nodes, oldPostion + 1, nodes, oldPostion, len); //TODO need this?
                    System.arraycopy(lineageCounts, oldPostion + 1, lineageCounts, oldPostion, len); //TODO need this?

                } else {
                    int len = oldPostion - newPosition;
                    // move everything from dest to src up one spot (opening up a gap at dest)
                    System.arraycopy(intervals, newPosition, intervals, newPosition + 1, len);
                    System.arraycopy(times, newPosition, times, newPosition + 1, len);
                    System.arraycopy(intervalTypes, newPosition, intervalTypes, newPosition + 1, len);
                    System.arraycopy(nodes, newPosition, nodes, newPosition + 1, len);//TODO need this?
                    System.arraycopy(lineageCounts, newPosition, lineageCounts, newPosition + 1, len);//TODO need this?
                }

                lineageCounts[newPosition] = lineageCount;
                intervals[newPosition] = interval;
                times[newPosition] = time;
                intervalTypes[newPosition] = type;
                nodes[newPosition] = node;

                // update node postions

                for (int i = Math.min(oldPostion, newPosition); i < Math.max(oldPostion, newPosition) + 1; i++) {
                    nodeOrder[nodes[i]] = i;
                }

                //update lineage counts
                // This current updated lineage counts separately from the changed intervals. Those could be done at the same time

                if (intervalTypes[newPosition] == IntervalType.COALESCENT) {
                    int difference = intervalTypes[newPosition - 1] == IntervalType.COALESCENT ? -1 : 1;
                    if (newPosition > oldPostion) {
                        for (int i = oldPostion; i < newPosition; i++) {
                             lineageCounts[i]++;
                        }
                        lineageCounts[newPosition] = lineageCounts[newPosition - 1] + difference;
                        // need to increment the lineage counts by 1
                    } else {
                        lineageCounts[newPosition] = lineageCounts[newPosition - 1] + difference;
                        for (int i = newPosition + 1; i < oldPostion + 1; i++) {
                            lineageCounts[i]--;
                        }
                    }
                } else { // was a sampling event
                    if (newPosition > oldPostion) {
                        int difference = intervalTypes[newPosition - 1] == IntervalType.COALESCENT ? -1 : 1;
                        for (int i = oldPostion; i < newPosition; i++) {
                            lineageCounts[i]--;
                        }
                        lineageCounts[newPosition] = lineageCounts[newPosition - 1] +difference;
                        // need to increment the lineage counts by 1
                    } else {

                        if (newPosition == 0) {
                            lineageCounts[newPosition] = 0;
                        } else {
                            int difference = intervalTypes[newPosition - 1] == IntervalType.COALESCENT ? -1 : 1;
                            lineageCounts[newPosition] = lineageCounts[newPosition-1] + difference;
                        }
                        for (int i = newPosition + 1; i < oldPostion + 1; i++) {
                            lineageCounts[i]++;
                        }
                    }
                }

                updateTimeAndIntervals(times[oldPostion], oldPostion);
            }


            updateTimeAndIntervals(newTime, newPosition);

        }

        public void copyEvents(Events source) {
            System.arraycopy(source.nodeOrder, 0, nodeOrder, 0, numberOfEvents);
            System.arraycopy(source.nodes, 0, nodes, 0, numberOfEvents);
            System.arraycopy(source.times, 0, times, 0, numberOfEvents);
            System.arraycopy(source.lineageCounts, 0, lineageCounts, 0, numberOfEvents);
            System.arraycopy(source.intervalTypes, 0, intervalTypes, 0, numberOfEvents);
            System.arraycopy(source.intervals, 0, intervals, 0, numberOfEvents);
        }

        /**
         * A pivate method that uses a binary search algorithm to find the first entry in the sorted event array
         * with a time greater than the proved value
         *
         * @param arr
         * @param value
         * @param startPos
         * @param endPos
         * @return
         */
        private int findFirstGreater(double[] arr, double value, int startPos, int endPos) {
            //https://stackoverflow.com/questions/6553970/find-the-first-element-in-a-sorted-array-that-is-greater-than-the-target
            //https://www.geeksforgeeks.org/first-strictly-greater-element-in-a-sorted-array-in-java/
            int low = startPos, high = endPos;
            int ans = -1;
            while (low <= high) {

                int mid = (low + high) / 2;
                if (arr[mid] <= value) {
                    /* This index, and everything below it, must not be the first element
                     * greater than what we're looking for because this element is no greater
                     * than the element.
                     */
                    low = mid + 1;
                } else {

                    ans = mid;
                    high = mid - 1;
                }

            }
            if (ans == -1) {
                // There is no greater element so the "first greater" is outside the array. This is handeled by
                // updateEventForNode
                ans = endPos + 1;
            }
            return ans;
        }

        public int getNode(int i) {
            return nodes[i];
        }

        public void updateTimeAndIntervals(double newTime, int i) {
            double newInterval = i == 0 ? -1 : newTime - times[i - 1];
            times[i] = newTime;
            intervals[i] = newInterval;

            if (i < numberOfEvents - 1) {
                double newNeighborInterval = times[i+1] - newTime;
                intervals[i+1] = newNeighborInterval;
            }
        }

        public void updateEventTime(double newTime, int i) {
            double newInterval = i == 0 ? -1 : newTime - times[i - 1];
            times[i] = newTime;
            intervals[i] = newInterval;
        }

        public int getLineageCount(int i) {
            return lineageCounts[i];
        }

        public double getInterval(int i) {
            return intervals[i];
        }

        public void setNodeOrder(int nodeNum, int position) {
            nodeOrder[nodeNum] = position;
        }


        private final int[] nodes;
        private final int[] nodeOrder;
        private final int[] lineageCounts;
        private final double[] intervals;
        private final double[] times;
        private final IntervalType[] intervalTypes;
        private final int numberOfEvents;


    }

    private final Events events;
    private final Events storedEvents;

    private List<Integer> updatedNodes;
    private List<Integer> storedUpdatedNodes;

    private boolean intervalsKnown;
    private boolean storedIntervalsKnown;

    private boolean onlyUpdateTimes;
    private boolean storedOnlyUpdateTimes;

    private final TreeModel tree;
    private boolean dirty;
    private int intervalCount = 0;


}

