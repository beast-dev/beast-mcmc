package dr.evomodel.bigfasttree;

import dr.evomodel.tree.TreeModel;

/**
 * @author Marc A. Suchard
 */
public class BestSignalsFromBigFastTreeIntervals extends BigFastTreeIntervals {

    private final Events originalEvents;

    public BestSignalsFromBigFastTreeIntervals(String name, TreeModel tree) {
        super(name, tree);
        originalEvents = new Events(tree.getNodeCount());
    }

    @Override
    public void calculateIntervals() {

        boolean changed = !intervalsKnown;

        if (!intervalsKnown) {
            originalEvents.copyEvents(events);
        }

        super.calculateIntervals();

        if (changed) {
            // Find first affected interval
            int firstEvent = compareEvents(originalEvents, events);
            if (firstEvent < events.size()) {
                fireModelChanged(new IntervalChangedEvent.FirstAffectedInterval(firstEvent - 1));
            }
        }
    }

    private int compareEvents(Events lhs, Events rhs) {

        int i = 0;
        for (; i < lhs.size(); ++i) {
            if (lhs.getNode(i) != rhs.getNode(i) || lhs.getInterval(i) != rhs.getInterval(i)) break;
        }

        return i;
    }
}
