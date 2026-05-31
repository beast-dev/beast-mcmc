package dr.inference.timeseries.representation;

/**
 * Read-only view of cached flat transition moments.
 *
 * <p>The arrays are owned by the transition representation. Callers may read them
 * until the underlying representation is marked dirty.</p>
 */
public final class TransitionMomentsView {

    private double[] transitionMatrix;
    private double[] transitionOffset;
    private double[] transitionCovariance;
    private int dimension;

    public void set(final double[] transitionMatrix,
                    final double[] transitionOffset,
                    final double[] transitionCovariance,
                    final int dimension) {
        this.transitionMatrix = transitionMatrix;
        this.transitionOffset = transitionOffset;
        this.transitionCovariance = transitionCovariance;
        this.dimension = dimension;
    }

    public double[] getTransitionMatrix() {
        return transitionMatrix;
    }

    public double[] getTransitionOffset() {
        return transitionOffset;
    }

    public double[] getTransitionCovariance() {
        return transitionCovariance;
    }

    public int getDimension() {
        return dimension;
    }
}
