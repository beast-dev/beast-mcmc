package dr.inference.hmc;

/**
 * Flat storage for a batch of per-parameter gradient vectors.
 */
public final class BatchGradient {

    private final int[] dimensions;
    private final int[] offsets;
    private final double[] values;

    public BatchGradient(final int[] dimensions) {
        if (dimensions == null || dimensions.length == 0) {
            throw new IllegalArgumentException("dimensions must contain at least one entry");
        }
        this.dimensions = dimensions.clone();
        this.offsets = new int[dimensions.length];
        int total = 0;
        for (int i = 0; i < dimensions.length; ++i) {
            if (dimensions[i] < 0) {
                throw new IllegalArgumentException("gradient dimension must be non-negative");
            }
            offsets[i] = total;
            total += dimensions[i];
        }
        this.values = new double[total];
    }

    public int size() {
        return dimensions.length;
    }

    public int dimension(final int index) {
        return dimensions[index];
    }

    public int offset(final int index) {
        return offsets[index];
    }

    public int totalDimension() {
        return values.length;
    }

    public double[] values() {
        return values;
    }

    public void setGradient(final int index, final double[] gradient) {
        checkGradient(index, gradient);
        System.arraycopy(gradient, 0, values, offsets[index], dimensions[index]);
    }

    public void copyGradientTo(final int index, final double[] destination, final int destinationOffset) {
        if (destination == null || destinationOffset < 0
                || destinationOffset + dimensions[index] > destination.length) {
            throw new IllegalArgumentException("gradient destination is too small");
        }
        System.arraycopy(values, offsets[index], destination, destinationOffset, dimensions[index]);
    }

    public void addGradient(final int index, final double[] increment) {
        checkGradient(index, increment);
        final int offset = offsets[index];
        final int dimension = dimensions[index];
        for (int i = 0; i < dimension; ++i) {
            values[offset + i] += increment[i];
        }
    }

    public void addBatch(final BatchGradient increment) {
        if (increment == null || increment.size() != size()) {
            throw new IllegalArgumentException("Gradient batch size mismatch");
        }
        for (int i = 0; i < size(); ++i) {
            if (increment.dimension(i) != dimensions[i]) {
                throw new IllegalArgumentException("Gradient dimension mismatch at batch index " + i);
            }
        }
        final double[] incrementValues = increment.values();
        for (int i = 0; i < values.length; ++i) {
            values[i] += incrementValues[i];
        }
    }

    private void checkGradient(final int index, final double[] gradient) {
        if (index < 0 || index >= dimensions.length) {
            throw new IllegalArgumentException("Gradient index out of bounds: " + index);
        }
        if (gradient == null || gradient.length != dimensions[index]) {
            throw new IllegalArgumentException("Gradient increment must have dimension " + dimensions[index]);
        }
    }
}
