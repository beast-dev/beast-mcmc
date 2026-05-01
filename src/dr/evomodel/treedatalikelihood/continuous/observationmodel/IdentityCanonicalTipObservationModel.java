package dr.evomodel.treedatalikelihood.continuous.observationmodel;

import dr.evomodel.treedatalikelihood.continuous.canonical.CanonicalTipObservation;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianState;

import java.util.Arrays;

public final class IdentityCanonicalTipObservationModel implements CanonicalTipObservationModel {

    private final int dimension;
    private final double[] values;
    private final boolean[] observed;
    private final int observedCount;

    public static IdentityCanonicalTipObservationModel fromObservation(
            final CanonicalTipObservation observation) {
        return new IdentityCanonicalTipObservationModel(
                observation.values,
                observation.observed,
                observation.observedCount);
    }

    public static IdentityCanonicalTipObservationModel missing(final int dimension) {
        return new IdentityCanonicalTipObservationModel(
                new double[dimension],
                new boolean[dimension],
                0);
    }

    private IdentityCanonicalTipObservationModel(final double[] values,
                                                final boolean[] observed,
                                                final int observedCount) {
        this.dimension = values.length;
        if (observed.length != dimension) {
            throw new IllegalArgumentException("observed mask dimension mismatch");
        }
        this.values = values.clone();
        this.observed = observed.clone();
        this.observedCount = observedCount;
    }

    @Override
    public int getLatentDimension() {
        return dimension;
    }

    @Override
    public int getObservationDimension() {
        return observedCount;
    }

    @Override
    public TipObservationMode getMode() {
        if (observedCount == 0) {
            return TipObservationMode.MISSING;
        }
        return observedCount == dimension
                ? TipObservationMode.EXACT_IDENTITY
                : TipObservationMode.PARTIAL_EXACT_IDENTITY;
    }

    @Override
    public boolean isEmpty() {
        return observedCount == 0;
    }

    public void copyTo(final CanonicalTipObservation observation) {
        if (observation.dim != dimension) {
            throw new IllegalArgumentException("Dimension mismatch");
        }
        System.arraycopy(values, 0, observation.values, 0, dimension);
        System.arraycopy(observed, 0, observation.observed, 0, dimension);
        observation.observedCount = observedCount;
    }

    public double[] getValues() {
        return values.clone();
    }

    public boolean[] getObserved() {
        return observed.clone();
    }

    public int getObservedCount() {
        return observedCount;
    }

    @Override
    public void fillChildCanonicalState(final CanonicalGaussianState out,
                                        final TipObservationModelWorkspace workspace) {
        throw new UnsupportedOperationException("Exact identity observations are handled by the collapsed path");
    }

    @Override
    public CanonicalTipObservationModel copy() {
        return new IdentityCanonicalTipObservationModel(values, observed, observedCount);
    }

    @Override
    public String toString() {
        return "IdentityCanonicalTipObservationModel{"
                + "mode=" + getMode()
                + ", values=" + Arrays.toString(values)
                + '}';
    }
}
