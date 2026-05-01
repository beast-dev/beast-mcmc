package dr.evomodel.treedatalikelihood.continuous.observationmodel;

import dr.evomodel.treedatalikelihood.preorder.NormalSufficientStatistics;
import org.ejml.data.DenseMatrix64F;

/**
 * Scratch-backed classifier for exact tip observations encoded as sufficient statistics.
 */
public final class SufficientStatisticsTipObservationPattern {

    private final int dimension;
    private final boolean[] observedMask;
    private final int[] observedIndices;
    private final int[] missingIndices;
    private final int[] reducedIndexByTrait;

    public SufficientStatisticsTipObservationPattern(final int dimension) {
        if (dimension < 1) {
            throw new IllegalArgumentException("dimension must be positive");
        }
        this.dimension = dimension;
        this.observedMask = new boolean[dimension];
        this.observedIndices = new int[dimension];
        this.missingIndices = new int[dimension];
        this.reducedIndexByTrait = new int[dimension];
    }

    public int classify(final NormalSufficientStatistics below) {
        final DenseMatrix64F precision = below.getRawPrecision();
        int observedCount = 0;
        for (int i = 0; i < dimension; ++i) {
            final double diagonal = precision.unsafe_get(i, i);
            if (Double.isInfinite(diagonal) && diagonal > 0.0) {
                observedMask[i] = true;
                observedCount++;
            } else if (diagonal == 0.0) {
                observedMask[i] = false;
            } else {
                return -1;
            }
            for (int j = 0; j < dimension; ++j) {
                if (i == j) {
                    continue;
                }
                if (precision.unsafe_get(i, j) != 0.0) {
                    return -1;
                }
            }
        }
        return observedCount;
    }

    public int collectPartition(final int observedCount) {
        int observedCursor = 0;
        int missingCursor = 0;
        for (int i = 0; i < dimension; ++i) {
            if (observedMask[i]) {
                observedIndices[observedCursor++] = i;
                reducedIndexByTrait[i] = -1;
            } else {
                missingIndices[missingCursor] = i;
                reducedIndexByTrait[i] = dimension + missingCursor;
                missingCursor++;
            }
        }
        if (observedCursor != observedCount) {
            throw new IllegalStateException("Observation partition count mismatch");
        }
        return missingCursor;
    }

    public boolean isObserved(final int trait) {
        return observedMask[trait];
    }

    public int observedIndex(final int offset) {
        return observedIndices[offset];
    }

    public int missingIndex(final int offset) {
        return missingIndices[offset];
    }

    public int reducedIndexByTrait(final int trait) {
        return reducedIndexByTrait[trait];
    }
}
