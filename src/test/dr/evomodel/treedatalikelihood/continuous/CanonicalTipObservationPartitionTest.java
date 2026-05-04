package test.dr.evomodel.treedatalikelihood.continuous;

import dr.evomodel.treedatalikelihood.continuous.observationmodel.CanonicalTipObservation;
import dr.evomodel.treedatalikelihood.continuous.observationmodel.IdentityCanonicalTipObservationModel;
import dr.evomodel.treedatalikelihood.continuous.observationmodel.TipObservationPartition;
import junit.framework.TestCase;

public final class CanonicalTipObservationPartitionTest extends TestCase {

    public void testPartitionCollectsObservedAndMissingIndices() {
        final CanonicalTipObservation observation = new CanonicalTipObservation(4);
        observation.setPartiallyObserved(
                new double[]{1.0, 2.0, 3.0, 4.0},
                new boolean[]{true, false, true, false});

        final TipObservationPartition partition = new TipObservationPartition(4);
        final int observedCount = partition.update(observation);

        assertEquals(2, observedCount);
        assertEquals(2, partition.getObservedCount());
        assertEquals(2, partition.getMissingCount());
        assertEquals(0, partition.observedIndex(0));
        assertEquals(2, partition.observedIndex(1));
        assertEquals(1, partition.missingIndex(0));
        assertEquals(3, partition.missingIndex(1));
        assertEquals(-1, partition.reducedIndexByTrait(0));
        assertEquals(4, partition.reducedIndexByTrait(1));
        assertEquals(-1, partition.reducedIndexByTrait(2));
        assertEquals(5, partition.reducedIndexByTrait(3));
    }

    public void testPartitionRejectsInconsistentObservedCount() {
        final CanonicalTipObservation observation = new CanonicalTipObservation(3);
        observation.setPartiallyObserved(
                new double[]{1.0, 2.0, 3.0},
                new boolean[]{true, false, true});
        observation.observedCount = 1;

        try {
            new TipObservationPartition(3).update(observation);
            fail("Expected inconsistent observedCount to be rejected");
        } catch (UnsupportedOperationException expected) {
            assertTrue(expected.getMessage().contains("observedCount"));
        }
    }

    public void testPartitionCollectsIdentityObservationModelIndices() {
        final CanonicalTipObservation observation = new CanonicalTipObservation(4);
        observation.setPartiallyObserved(
                new double[]{1.0, 2.0, 3.0, 4.0},
                new boolean[]{false, true, false, true});
        final IdentityCanonicalTipObservationModel model =
                IdentityCanonicalTipObservationModel.fromObservation(observation);

        final TipObservationPartition partition = new TipObservationPartition(4);
        final int observedCount = partition.update(model);

        assertEquals(2, observedCount);
        assertEquals(1, partition.observedIndex(0));
        assertEquals(3, partition.observedIndex(1));
        assertEquals(0, partition.missingIndex(0));
        assertEquals(2, partition.missingIndex(1));
        assertEquals(4, partition.reducedIndexByTrait(0));
        assertEquals(-1, partition.reducedIndexByTrait(1));
        assertEquals(5, partition.reducedIndexByTrait(2));
        assertEquals(-1, partition.reducedIndexByTrait(3));
    }
}
