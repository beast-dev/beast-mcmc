package dr.evolution.coalescent;

public interface IntervalEventList extends IntervalList {
    /**
	 * Returns the number of events in the interval list.
	 * should be 1 more than interval count.
	 */
	int getEventCount();

	/**
	 * Returns the time of the ith event.
	 */
	double getEventTime(int i);
}
