package dr.inference.timeseries.core;

/**
 * Time grid abstraction for discrete or irregularly spaced observation times.
 */
public interface TimeGrid {

    int getTimeCount();

    double getTime(int index);

    double getDelta(int fromIndex, int toIndex);

    boolean isRegular();
}
