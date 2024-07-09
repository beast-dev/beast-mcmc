package dr.evomodel.bigfasttree;

import dr.evomodel.tree.TreeModel;

/**
 * @author Marc A. Suchard
 * @author Guy Baele
 */
public class BestSignalsFromBigFastTreeIntervals extends BigFastTreeIntervals {

    private final Events originalEvents;

    public BestSignalsFromBigFastTreeIntervals(TreeModel tree) {
        super("signalsFromBigFastIntervals", tree);
        originalEvents = new Events(tree.getNodeCount());
    }

    public BestSignalsFromBigFastTreeIntervals(String name, TreeModel tree) {
        super(name, tree);
        originalEvents = new Events(tree.getNodeCount());
    }

    @Override
    public void calculateIntervals() {

        boolean changed = !intervalsKnown;

        if (!intervalsKnown && (originalEvents != null)) {
            originalEvents.copyEvents(events);
        }

        super.calculateIntervals();

        if (changed && (originalEvents != null)) {
            // Find first affected interval
            /*int firstEvent = findFirstEvent(originalEvents, events);
            if (firstEvent < events.size()) {
                fireModelChanged(new IntervalChangedEvent.FirstAffectedInterval(firstEvent - 1));
            }*/

            //using the code below to return range of intervals that are changed
            //find first and last interval affected
            int[] changedEvents = findChangedEvents(originalEvents, events);
            if (changedEvents[0] < events.size() && changedEvents[1] < events.size()) {
                fireModelChanged(new IntervalChangedEvent.AffectedIntervals(changedEvents));
            }
        }
    }

    private int findFirstEvent(Events lhs, Events rhs) {
        int i = 0;
        for (; i < lhs.size(); ++i) {
            if (lhs.getNode(i) != rhs.getNode(i) || lhs.getInterval(i) != rhs.getInterval(i)) break;
        }
        return i;
    }

    private int[] findChangedEvents(Events lhs, Events rhs) {
        int[] changes = new int[2];
        changes[0] = -1;
        changes[1] = -1;
        for (int i = 0; i < lhs.size(); ++i) {
            if (lhs.getNode(i) != rhs.getNode(i) || lhs.getInterval(i) != rhs.getInterval(i)) {
                if (changes[0] == -1) {
                    changes[0] = i-1;
                } else {
                    changes[1] = i-1;
                }
            }
        }
        return changes;
    }
}
