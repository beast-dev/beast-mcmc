package dr.evomodel.treedatalikelihood.continuous.observationmodel;

import dr.evomodel.treedatalikelihood.continuous.canonical.CanonicalTipObservation;

public final class TipObservationPartition {

    public final int[] observedIndices;
    public final int[] missingIndices;
    public final int[] reducedIndexByTrait;
    private final int dimension;
    private int observedCount;
    private int missingCount;

    public TipObservationPartition(final int dimension) {
        if (dimension < 1) {
            throw new IllegalArgumentException("dimension must be positive");
        }
        this.dimension = dimension;
        this.observedIndices = new int[dimension];
        this.missingIndices = new int[dimension];
        this.reducedIndexByTrait = new int[dimension];
    }

    public int update(final CanonicalTipObservation observation) {
        if (observation.dim != dimension) {
            throw new IllegalArgumentException("Observation dimension mismatch");
        }

        observedCount = 0;
        missingCount = 0;
        for (int i = 0; i < dimension; ++i) {
            if (observation.observed[i]) {
                observedIndices[observedCount++] = i;
                reducedIndexByTrait[i] = -1;
            } else {
                missingIndices[missingCount++] = i;
                reducedIndexByTrait[i] = dimension + missingCount - 1;
            }
        }

        if (observedCount != observation.observedCount || observedCount + missingCount != dimension) {
            throw new UnsupportedOperationException(
                    "Canonical tip observation partition is inconsistent with observedCount.");
        }
        return observedCount;
    }

    public int getDimension() {
        return dimension;
    }

    public int getObservedCount() {
        return observedCount;
    }

    public int getMissingCount() {
        return missingCount;
    }

    public int observedIndex(final int index) {
        return observedIndices[index];
    }

    public int missingIndex(final int index) {
        return missingIndices[index];
    }

    public int reducedIndexByTrait(final int trait) {
        return reducedIndexByTrait[trait];
    }

    public boolean isObserved(final CanonicalTipObservation observation, final int trait) {
        return observation.observed[trait];
    }
}
