package dr.inference.timeseries.core;

/**
 * Immutable time grid backed by explicit, strictly increasing observation times.
 */
public final class IrregularTimeGrid implements TimeGrid {

    private final double[] times;

    public IrregularTimeGrid(final double[] times) {
        if (times == null || times.length < 1) {
            throw new IllegalArgumentException("times must contain at least one entry");
        }
        this.times = times.clone();
        for (int i = 1; i < this.times.length; ++i) {
            if (!(this.times[i] > this.times[i - 1])) {
                throw new IllegalArgumentException("times must be strictly increasing");
            }
        }
    }

    @Override
    public int getTimeCount() {
        return times.length;
    }

    @Override
    public double getTime(final int index) {
        checkIndex(index);
        return times[index];
    }

    @Override
    public double getDelta(final int fromIndex, final int toIndex) {
        checkIndex(fromIndex);
        checkIndex(toIndex);
        if (toIndex < fromIndex) {
            throw new IllegalArgumentException("toIndex must be >= fromIndex");
        }
        return times[toIndex] - times[fromIndex];
    }

    @Override
    public boolean isRegular() {
        return false;
    }

    public double[] getTimes() {
        return times.clone();
    }

    private void checkIndex(final int index) {
        if (index < 0 || index >= times.length) {
            throw new IndexOutOfBoundsException("Invalid time index: " + index);
        }
    }
}
