package dr.evomodel.treedatalikelihood.continuous;

import dr.evomodel.treedatalikelihood.preorder.NormalSufficientStatistics;
import org.ejml.data.DenseMatrix64F;

/**
 * Scratch-backed classifier for exact canonical tip observations.
 */
final class OUCanonicalObservationPattern {

    private final int dimension;
    private final boolean[] observedMask;
    private final int[] observedIndices;
    private final int[] missingIndices;
    private final int[] reducedIndexByTrait;

    OUCanonicalObservationPattern(final int dimension) {
        this.dimension = dimension;
        this.observedMask = new boolean[dimension];
        this.observedIndices = new int[dimension];
        this.missingIndices = new int[dimension];
        this.reducedIndexByTrait = new int[dimension];
    }

    int classify(final NormalSufficientStatistics below) {
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

    int collectPartition(final int observedCount) {
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

    boolean isObserved(final int trait) {
        return observedMask[trait];
    }

    int observedIndex(final int offset) {
        return observedIndices[offset];
    }

    int missingIndex(final int offset) {
        return missingIndices[offset];
    }

    int reducedIndexByTrait(final int trait) {
        return reducedIndexByTrait[trait];
    }
}
