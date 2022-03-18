package dr.evomodel.bigfasttree;

/**
 * @author Marc A. Suchard
 */
public interface IntervalChangedEvent {

    int getInterval();

    class FirstAffectedInterval implements IntervalChangedEvent {

        private final int interval;

        public FirstAffectedInterval(int interval) {
            this.interval = interval;
        }

        @Override
        public int getInterval() {
            return interval;
        }
    }
}
