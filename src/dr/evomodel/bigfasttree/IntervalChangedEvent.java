package dr.evomodel.bigfasttree;

/**
 * @author Marc A. Suchard
 */
public interface IntervalChangedEvent {

    int getInterval();

    class FirstAffectedInterval implements IntervalChangedEvent {

        //TODO expand this to first and last intervals, to minimize redoing the number of matrix exponentiations?
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
