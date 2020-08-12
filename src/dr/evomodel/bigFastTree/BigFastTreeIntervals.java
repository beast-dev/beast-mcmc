package dr.evomodel.bigFastTree;

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


    public void makeDirty(){
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
        return events.getInterval(i+1);
    }

    @Override
    public double getIntervalTime(int i){
        if (!intervalsKnown){
            calculateIntervals();
        }
        return events.getTime(i);
    }

    @Override
    public int getLineageCount(int i) {
        if (!intervalsKnown) {
            calculateIntervals();
        }
        return events.getLineageCount(i+1);
    }

    @Override
    public int getCoalescentEvents(int i) {
        if (!intervalsKnown) {
            calculateIntervals();
        }
        if (i < intervalCount - 1) {
            return events.getLineageCount(i+1) - events.getLineageCount(i+2);
        } else {
            return events.getLineageCount(i+1) - 1;
        }
    }

    @Override
    public IntervalType getIntervalType(int i) {
        if (!intervalsKnown) {
            calculateIntervals();
        }
        return events.getType(i+1);
    }

    @Override
    public double getTotalDuration() {
        if (!intervalsKnown) {
            calculateIntervals();
        }
        return events.getTime(events.getSize()-1);
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
        //If dirty we rebuild the evens and sort them using parrelle sort
        if (dirty) {
            // Resort nodes by heights
            // will this update the tree nodes?

            NodeRef[] nodes = new NodeRef[tree.getNodeCount()];
            System.arraycopy(tree.getNodes(),0,nodes,0,tree.getNodeCount());
            Arrays.parallelSort(nodes, (a, b) -> Double.compare(tree.getNodeHeight(a), tree.getNodeHeight(b)));

            intervalCount = nodes.length - 1;

            double lastTime = tree.getNodeHeight(nodes[0]);
            if (!tree.isExternal(nodes[0])) {
                throw new IllegalArgumentException("The first event is not a sample event");
            }
            events.setEvent(new Event(lastTime,IntervalType.SAMPLE,nodes[0].getNumber(),-1,0),0);

            int lineages = 1;
            for (int i = 1; i < nodes.length; i++) {
                NodeRef node = nodes[i];
                double time = tree.getNodeHeight(node);
                double interval = time - lastTime;
                IntervalType type = tree.isExternal(node) ? IntervalType.SAMPLE : IntervalType.COALESCENT;
                int lineageCount = lineages;
                events.setEvent(new Event(time,type,node.getNumber(),interval,lineageCount),i);
                events.setNodeOrder(node.getNumber(), i);


                if (type == IntervalType.SAMPLE) {
                    lineages++;
                } else {
                    lineages--;
                }
                lastTime = time;
            }
            intervalsKnown = true;

        }else if (onlyUpdateTimes) {
            for (int i = 0; i < events.getSize(); i++) {
                double newTime = tree.getNodeHeight(tree.getNode(events.getNode(i)));
                events.updateTime(newTime, i);
            }
            onlyUpdateTimes=false;
        }else {
            for(int node : updatedNodes){
                events.updateForChangedNode(node, tree.getNodeHeight(tree.getNode(node)));
            }
        }

        intervalsKnown = true;
        dirty=false;
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
                        intervalsKnown= false;
                    }

                }else if (treeChangedEvent.isTreeChanged()) {
                    if(!treeChangedEvent.isNodeOrderChanged()){
                        onlyUpdateTimes=true;
                    }else {
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
        Collections.copy(storedUpdatedNodes,updatedNodes);
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

        Event(double time, IntervalType type, int node,double interval, int lineageCount) {
            this.time = time;
            this.type = type;
            this.node = node;
            this.interval = interval;
            this.lineageCount = lineageCount;
//            this.info = info;
        }

         /* The type of event
         */
        final IntervalType  type;

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

        public Events(int n) {
            events = new Event[n];
            nodeOrder = new int[n];
        }


        public int getSize(){
            return events.length;
        }
        public void sortEvents() {
            Arrays.parallelSort(events);
            for (int i = 0; i < events.length; i++) {
                nodeOrder[events[i].node] = i;
            }
        }


        public int getNodePosition(int nodeNum) {
            return nodeOrder[nodeNum];
        }

        public IntervalType getType(int i) {
            return events[i].type;
        }

        public double getTime(int i) {
            return events[i].time;
        }

        public int size() {
            return events.length;
        }

        public void setEvent(Event e, int i) {
            events[i] = e;
        }

        public void updateForChangedNode(int node, double newTime) {
            int oldPostion = nodeOrder[node];
            double oldHieght = events[oldPostion].time;

            int newPosition=oldPostion;
                if (newTime > oldHieght) {
                    int firstHeigher = findFirstGreater(events, newTime, oldPostion, events.length-1);
                    newPosition = firstHeigher-1;

                }else if(newTime <oldHieght){
                    newPosition = findFirstGreater(events, newTime, 0, oldPostion-1);
                }


            if(newPosition!=oldPostion) {
                // Borrowed from SortableArray
                    Event event = events[oldPostion];
                    if (newPosition > oldPostion) {
                        int len = newPosition - oldPostion;
                        // move everything from the src to the dest down one spot (opening up a gap at dest)
                        System.arraycopy(events, oldPostion + 1, events, oldPostion, len);
                    } else {
                        int len = oldPostion - newPosition;
                        // move everything from dest to src up one spot (opening up a gap at dest)
                        System.arraycopy(events, newPosition, events, newPosition + 1, len);
                    }
                    events[newPosition] = event;

                    // update node postions

                    for (int i = Math.min(oldPostion, newPosition); i < Math.max(oldPostion, newPosition)+1; i++) {
                        nodeOrder[events[i].node] = i;
                    }

                //update lineage counts
                // This current updated lineage counts separately from the changed intervals. Those could be done at the same time

                if(events[newPosition].type==IntervalType.COALESCENT){
                    int difference = events[newPosition-1].type==IntervalType.COALESCENT?-1: 1;
                    if(newPosition>oldPostion){
                        for (int i = oldPostion; i < newPosition ; i++) {
                            Event e = events[i];
                            events[i] = new Event(e.time, e.type, e.node, e.interval, e.lineageCount + 1);
                        }
                        Event movedEvent = events[newPosition];
                        events[newPosition] = new Event(movedEvent.time, movedEvent.type, movedEvent.node, movedEvent.interval, events[newPosition-1].lineageCount + difference);
                        // need to increment the lineage counts by 1
                    }else{
                        Event movedEvent = events[newPosition];
                        events[newPosition] = new Event(movedEvent.time, movedEvent.type, movedEvent.node, movedEvent.interval, events[newPosition-1].lineageCount + difference);
                        for (int i = newPosition+1; i < oldPostion+1 ; i++) {
                            Event e = events[i];
                            events[i] = new Event(e.time, e.type, e.node, e.interval, e.lineageCount - 1);
                        }
                    }
                }else{ // was a sampling event
                    if(newPosition>oldPostion){
                        int difference = events[newPosition-1].type==IntervalType.COALESCENT?-1: 1;
                        for (int i = oldPostion; i < newPosition ; i++) {
                            Event e = events[i];
                            events[i] = new Event(e.time, e.type, e.node, e.interval, e.lineageCount - 1);
                        }
                        Event movedEvent = events[newPosition];
                        events[newPosition] = new Event(movedEvent.time, movedEvent.type, movedEvent.node, movedEvent.interval, events[newPosition-1].lineageCount +difference);
                        // need to increment the lineage counts by 1
                    }else{
                        Event movedEvent = events[newPosition];
                        if(newPosition==0){
                            events[newPosition] = new Event(movedEvent.time, movedEvent.type, movedEvent.node, -1, 0);
                        }else{
                            int difference = events[newPosition-1].type==IntervalType.COALESCENT?-1: 1;
                            events[newPosition] = new Event(movedEvent.time, movedEvent.type, movedEvent.node, movedEvent.interval, events[newPosition-1].lineageCount + difference);
                        }
                        for (int i = newPosition+1; i < oldPostion +1; i++) {
                            Event e = events[i];
                            events[i] = new Event(e.time, e.type, e.node, e.interval, e.lineageCount + 1);
                        }
                    }
                }
                updateTime(events[oldPostion].time, oldPostion);
            }

            updateTime(newTime,newPosition);

        }

        public void copyEvents(Events source) {
            System.arraycopy(source.events,0,events,0,source.events.length);
            System.arraycopy(source.nodeOrder,0,nodeOrder,0,source.nodeOrder.length);
        }

        /**
         * A pivate method that uses a binary search algorithm to find the first entry in the sorted event array
         * with a time greater than the proved value
         * @param arr
         * @param value
         * @param startPos
         * @param endPos
         * @return
         */
        private int findFirstGreater(Event[] arr, double value, int startPos, int endPos) {
            //https://stackoverflow.com/questions/6553970/find-the-first-element-in-a-sorted-array-that-is-greater-than-the-target
            //https://www.geeksforgeeks.org/first-strictly-greater-element-in-a-sorted-array-in-java/
            int low = startPos, high = endPos;
            int ans = -1;
            while (low <= high) {

                int mid = (low + high) / 2;
                if (arr[mid].time <= value) {
                    /* This index, and everything below it, must not be the first element
                     * greater than what we're looking for because this element is no greater
                     * than the element.
                     */
                    low = mid + 1;
                }
                else {

                    ans = mid;
                    high = mid - 1;
                }

            }
            if(ans==-1){
                // There is no greater element so the "first greater" is outside the array. This is handeled by
                // updateEventForNode
                ans=endPos+1;
            }
            return ans;
        }

        private final Event[] events;
        private final int[] nodeOrder;

        public int getNode(int i) {
            return events[i].node;
        }

        public void updateTime(double newTime, int i) {
            Event oldEvent = events[i];
            double newInterval = i==0? -1 : newTime-events[i-1].time;
            events[i] = new Event(newTime, oldEvent.type, oldEvent.node, newInterval, oldEvent.lineageCount);

            if(i<events.length-1) {
                Event oldNeighbor = events[i + 1];
                double newNeighborInterval = oldNeighbor.time - newTime;
                events[i + 1] = new Event(oldNeighbor.time, oldNeighbor.type, oldNeighbor.node, newNeighborInterval, oldNeighbor.lineageCount);
            }
        }

        public int getLineageCount(int i) {
            return events[i].lineageCount;
        }

        public double getInterval(int i) {
            return events[i].interval;
        }

        public void setNodeOrder(int nodeNum, int position) {
            nodeOrder[nodeNum]=position;
        }
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

