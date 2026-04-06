package dr.inference.model;

import java.util.Objects;
/*
 * IndexedParameter.java
 *
 * Deterministic  Parameter:
 *
 *   out[b] = values[ round(indices[b]) - indexOffset ]
 *
 * Typical use:
 *   - values  (dimension K, real in [0,1])
 *   - indices (dimension B, integer-like in {1..K} or {0..K-1})
 *
 * This parameter is IMMUTABLE: operators should act on "values" and "indices",
 * never on this derived parameter.
 */

/**
 * @author Filippo Monti
 */

public final class IndexedParameter extends Parameter.Abstract implements VariableListener {

    private final Parameter values;     // lookup table
    private final Parameter indices;    // index vector (integer-like doubles)
    private final int indexOffset;      // 0 if indices are 0-based; 1 if indices are 1-based
    private final Bounds<Double> bounds;

    /**
     * @param id          object id
     * @param values      lookup parameter (dimension K)
     * @param indices     index parameter (dimension B), integer-like
     * @param indexOffset 0 for 0-based indices, 1 for 1-based indices
     */
    public IndexedParameter(final String id,
                            final Parameter values,
                            final Parameter indices,
                            final int indexOffset) {
        super(id);

        this.values = Objects.requireNonNull(values, "values");
        this.indices = Objects.requireNonNull(indices, "indices");
        this.indexOffset = indexOffset;

        // Propagate changes from upstream parameters.
        this.values.addParameterListener(this);
        this.indices.addParameterListener(this);

        // Bounds: replicate a conservative bound from "values" onto each derived dimension.
        this.bounds = makeReplicatedBounds(values, indices.getDimension());
    }

    public Parameter getValuesParameter() {
        return values;
    }

    public Parameter getIndicesParameter() {
        return indices;
    }

    public int getIndexOffset() {
        return indexOffset;
    }

    @Override
    public int getDimension() {
        return indices.getDimension();
    }

    @Override
    public double getParameterValue(final int dim) {
        final double raw = indices.getParameterValue(dim);
        final int idx = toIndexStrict(raw) - indexOffset;

        if (idx < 0 || idx >= values.getDimension()) {
            throw new IllegalArgumentException(
                    "IndexedParameter '" + getId() + "': index out of range at dim=" + dim +
                            " (raw=" + raw + ", indexOffset=" + indexOffset +
                            ", resolved idx=" + idx + ", valuesDim=" + values.getDimension() + ")"
            );
        }
        return values.getParameterValue(idx);
    }

    public int getIndexValue(final int dim) {
        final double raw = indices.getParameterValue(dim);
        return toIndexStrict(raw);
    }

    @Override
    public String getParameterName() {
        return getId();
    }

    @Override
    public Bounds<Double> getBounds() {
        return bounds;
    }

    @Override
    public void addBounds(final Bounds<Double> bounds) {
        // Deterministic parameter: bounds should be specified on "values" instead.
        throw new UnsupportedOperationException(
                "IndexedParameter is deterministic; add bounds to the underlying 'values' parameter."
        );
    }

    @Override
    public void setDimension(final int dim) {
        throw new UnsupportedOperationException("IndexedParameter is deterministic; dimension is indices.getDimension().");
    }

    @Override
    public void addDimension(final int index, final double value) {
        throw new UnsupportedOperationException("IndexedParameter is deterministic; update the 'indices' parameter.");
    }

    @Override
    public double removeDimension(final int index) {
        throw new UnsupportedOperationException("IndexedParameter is deterministic; update the 'indices' parameter.");
    }

    // -----------------------
    // Mutators: DISALLOWED
    // -----------------------

    @Override
    public void setParameterValue(final int dim, final double value) {
        throw new UnsupportedOperationException("IndexedParameter is deterministic; update 'values'/'indices' instead.");
    }

    @Override
    public void setParameterValueQuietly(final int dim, final double value) {
        throw new UnsupportedOperationException("IndexedParameter is deterministic; update 'values'/'indices' instead.");
    }

    @Override
    public void setParameterValueNotifyChangedAll(final int dim, final double value) {
        throw new UnsupportedOperationException("IndexedParameter is deterministic; update 'values'/'indices' instead.");
    }

    @Override
    public boolean isImmutable() {
        return true;
    }

    // -----------------------
    // State handling: no-op (derived)
    // -----------------------

    @Override
    protected void storeValues() { /* no-op */ }

    @Override
    protected void restoreValues() { /* no-op */ }

    @Override
    protected void acceptValues() { /* no-op */ }

    @Override
    protected void adoptValues(final Parameter source) { /* no-op */ }

    // -----------------------
    // Upstream change propagation
    // -----------------------

    @Override
    public void variableChangedEvent(final Variable variable, final int index, final Parameter.ChangeType type) {
        // Any upstream change can affect any derived dimension => fire "all changed".
        fireParameterChangedEvent(-1, Parameter.ChangeType.ALL_VALUES_CHANGED);
    }

    // -----------------------
    // Helpers
    // -----------------------

    private static int toIndexStrict(final double x) {
        // enforce integer-like indices.
        final long r = Math.round(x);
        if (Math.abs(x - r) > 1e-9) {
            throw new IllegalArgumentException("IndexedParameter: index must be integer-like but found " + x);
        }
        if (r < Integer.MIN_VALUE || r > Integer.MAX_VALUE) {
            throw new IllegalArgumentException("IndexedParameter: index out of int range: " + x);
        }
        return (int) r;
    }

    private static Bounds<Double> makeReplicatedBounds(final Parameter values, final int outDim) {
        double lower = Double.NEGATIVE_INFINITY;
        double upper = Double.POSITIVE_INFINITY;

        try {
            final Bounds<Double> vb = values.getBounds();
            // Use bounds at dimension 0 as representative.
            lower = vb.getLowerLimit(0);
            upper = vb.getUpperLimit(0);
        } catch (Throwable ignored) {
            // leave infinities
        }

        final double lo = lower;
        final double hi = upper;

        return new Bounds<Double>() {
            @Override
            public Double getUpperLimit(int i) { return hi; }

            @Override
            public Double getLowerLimit(int i) { return lo; }

            @Override
            public int getBoundsDimension() { return outDim; }

        };
    }
}