/*
 * StructuredCoalescentScheduleWalker.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
 */

package dr.evomodel.coalescent;

/**
 * Shared biological event loop for structured-coalescent engines.
 *
 * The visitor sees the active-lineage set before each event is applied. The
 * walker then updates the active set itself, so MASCOT and BASTA cannot drift
 * in sample/coalescent ordering semantics.
 */
public final class StructuredCoalescentScheduleWalker {

    public interface Visitor {
        void initialSample(int sampleNode, StructuredCoalescentActiveLineages activeLineages);

        void interval(int interval, double intervalLength, int eventType,
                      StructuredCoalescentActiveLineages activeLineages);

        void sampleEvent(int interval, double intervalLength, int sampleNode,
                         StructuredCoalescentActiveLineages activeLineages);

        void coalescentEvent(int interval, double intervalLength, int leftChild, int rightChild, int parent,
                             StructuredCoalescentActiveLineages activeLineages);

        void afterEvent(int interval, int eventType, StructuredCoalescentActiveLineages activeLineages);

        void finish(StructuredCoalescentActiveLineages activeLineages);
    }

    public static class Adapter implements Visitor {
        @Override
        public void initialSample(int sampleNode, StructuredCoalescentActiveLineages activeLineages) {
        }

        @Override
        public void interval(int interval, double intervalLength, int eventType,
                             StructuredCoalescentActiveLineages activeLineages) {
        }

        @Override
        public void sampleEvent(int interval, double intervalLength, int sampleNode,
                                StructuredCoalescentActiveLineages activeLineages) {
        }

        @Override
        public void coalescentEvent(int interval, double intervalLength, int leftChild, int rightChild, int parent,
                                    StructuredCoalescentActiveLineages activeLineages) {
        }

        @Override
        public void afterEvent(int interval, int eventType, StructuredCoalescentActiveLineages activeLineages) {
        }

        @Override
        public void finish(StructuredCoalescentActiveLineages activeLineages) {
        }
    }

    private StructuredCoalescentScheduleWalker() {
    }

    public static void validate(StructuredCoalescentSchedule schedule) {
        walk(schedule, new Adapter());
    }

    public static void walk(StructuredCoalescentSchedule schedule, Visitor visitor) {
        if (schedule == null) {
            throw new IllegalArgumentException("schedule must not be null");
        }
        if (visitor == null) {
            throw new IllegalArgumentException("visitor must not be null");
        }

        StructuredCoalescentActiveLineages activeLineages =
                new StructuredCoalescentActiveLineages(schedule.nodeCount);
        activeLineages.reset(schedule.initialSampleNode);
        visitor.initialSample(schedule.initialSampleNode, activeLineages);

        for (int interval = 0; interval < schedule.getIntervalCount(); interval++) {
            double intervalLength = schedule.intervalLengths[interval];
            int eventType = schedule.eventTypes[interval];
            visitor.interval(interval, intervalLength, eventType, activeLineages);

            if (eventType == StructuredCoalescentSchedule.SAMPLE) {
                int sampleNode = schedule.sampleNodes[interval];
                visitor.sampleEvent(interval, intervalLength, sampleNode, activeLineages);
                activeLineages.addSample(sampleNode);
            } else if (eventType == StructuredCoalescentSchedule.COALESCENT) {
                int leftChild = schedule.coalescentChild1[interval];
                int rightChild = schedule.coalescentChild2[interval];
                int parent = schedule.coalescentParent[interval];
                visitor.coalescentEvent(interval, intervalLength, leftChild, rightChild, parent, activeLineages);
                activeLineages.replaceChildrenWithParent(leftChild, rightChild, parent);
            } else {
                throw new IllegalArgumentException("unknown event type at " + interval + ": " + eventType);
            }

            visitor.afterEvent(interval, eventType, activeLineages);
        }

        activeLineages.requireSingleActiveLineage();
        visitor.finish(activeLineages);
    }
}
