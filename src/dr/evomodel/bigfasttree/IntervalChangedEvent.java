package dr.evomodel.bigfasttree;

/**
 * @author Marc A. Suchard
 * @author Guy Baele
 */
public interface IntervalChangedEvent {

    int getFirstInterval();
    int getLastInterval();

    class FirstAffectedInterval implements IntervalChangedEvent {

        private final int firstInterval;

        public FirstAffectedInterval(int interval) {
            this.firstInterval = interval;
        }

        @Override
        public int getFirstInterval() {
            return firstInterval;
        }

        @Override
        public int getLastInterval() {
            throw new RuntimeException("");
        }

    }

    class AffectedIntervals implements IntervalChangedEvent {

        private final int firstInterval;
        private final int lastInterval;

        public AffectedIntervals(int firstInterval, int lastInterval) {
            this.firstInterval = firstInterval;
            this.lastInterval = lastInterval;
        }

        public AffectedIntervals(int[] intervals) {
            this.firstInterval = intervals[0];
            this.lastInterval = intervals[1];
        }

        @Override
        public int getFirstInterval() {
            return firstInterval;
        }

        @Override
        public int getLastInterval() {
            return lastInterval;
        }

    }

}
