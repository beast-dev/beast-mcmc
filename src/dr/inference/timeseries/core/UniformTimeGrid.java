package dr.inference.timeseries.core;

/**
 * Simple immutable uniform time grid.
 */
public final class UniformTimeGrid implements TimeGrid {

    private final int timeCount;
    private final double startTime;
    private final double timeStep;

    public UniformTimeGrid(final int timeCount,
                           final double startTime,
                           final double timeStep) {
        if (timeCount < 1) {
            throw new IllegalArgumentException("timeCount must be at least 1");
        }
        if (timeStep <= 0.0) {
            throw new IllegalArgumentException("timeStep must be > 0");
        }
        this.timeCount = timeCount;
        this.startTime = startTime;
        this.timeStep = timeStep;
    }

    @Override
    public int getTimeCount() {
        return timeCount;
    }

    @Override
    public double getTime(final int index) {
        checkIndex(index);
        return startTime + index * timeStep;
    }

    @Override
    public double getDelta(final int fromIndex, final int toIndex) {
        checkIndex(fromIndex);
        checkIndex(toIndex);
        if (toIndex < fromIndex) {
            throw new IllegalArgumentException("toIndex must be >= fromIndex");
        }
        return (toIndex - fromIndex) * timeStep;
    }

    @Override
    public boolean isRegular() {
        return true;
    }

    public double getStartTime() {
        return startTime;
    }

    public double getTimeStep() {
        return timeStep;
    }

    private void checkIndex(final int index) {
        if (index < 0 || index >= timeCount) {
            throw new IndexOutOfBoundsException("Invalid time index: " + index);
        }
    }
}
