/*
 * StructuredCoalescentTipDataTest.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
 */

package test.dr.evomodel.coalescent;

import dr.evomodel.coalescent.StructuredCoalescentSchedule;
import dr.evomodel.coalescent.StructuredCoalescentTipData;
import junit.framework.TestCase;

import java.util.Arrays;

public class StructuredCoalescentTipDataTest extends TestCase {

    public void testNormalizesPartialsAndWritesObservedState() {
        double[][] partials = new double[3][];
        partials[0] = new double[]{2.0, 1.0, 0.0};
        int[] observedStates = filledObservedStates(3);
        observedStates[1] = 2;

        StructuredCoalescentTipData tipData = new StructuredCoalescentTipData(3, 3, partials, observedStates);
        tipData.validateCompatibleWith(twoTipSchedule());

        assertTrue(tipData.hasTipPartials(0));
        assertFalse(tipData.hasTipPartials(1));
        assertFalse(tipData.hasObservedState(0));
        assertTrue(tipData.hasObservedState(1));
        assertDoubleArray(new double[]{2.0 / 3.0, 1.0 / 3.0, 0.0}, tipData.getTipPartialsCopy(0), 0.0);
        assertDoubleArray(new double[]{0.0, 0.0, 1.0}, tipData.getTipPartialsCopy(1), 0.0);

        double[] out = new double[]{-1.0, -1.0, -1.0, -1.0, -1.0};
        tipData.writeTipPartials(1, out, 1);
        assertDoubleArray(new double[]{-1.0, 0.0, 0.0, 1.0, -1.0}, out, 0.0);
    }

    public void testCopiesInputArrays() {
        double[][] partials = new double[3][];
        partials[0] = new double[]{1.0, 1.0};
        int[] observedStates = filledObservedStates(3);
        observedStates[1] = 1;

        StructuredCoalescentTipData tipData = new StructuredCoalescentTipData(3, 2, partials, observedStates);
        partials[0][0] = 100.0;
        observedStates[1] = 0;

        assertDoubleArray(new double[]{0.5, 0.5}, tipData.getTipPartialsCopy(0), 0.0);
        assertDoubleArray(new double[]{0.0, 1.0}, tipData.getTipPartialsCopy(1), 0.0);
    }

    public void testValidateRejectsMissingSamplePayload() {
        double[][] partials = new double[3][];
        partials[1] = new double[]{1.0, 0.0};
        StructuredCoalescentTipData tipData = StructuredCoalescentTipData.fromPartials(3, 2, partials);

        try {
            tipData.validateCompatibleWith(twoTipSchedule());
            fail("expected missing initial sample payload to be rejected");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage().contains("missing tip state or partials"));
        }
    }

    public void testRejectsInvalidPartials() {
        try {
            StructuredCoalescentTipData.fromPartials(2, 2, new double[][]{new double[]{0.0, 0.0}});
            fail("expected zero-sum partials to be rejected");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage().contains("positive finite sum"));
        }

        try {
            StructuredCoalescentTipData.fromPartials(2, 2, new double[][]{new double[]{1.0, -0.1}});
            fail("expected negative partials to be rejected");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage().contains("invalid tip partial"));
        }
    }

    private static StructuredCoalescentSchedule twoTipSchedule() {
        return new StructuredCoalescentSchedule(
                0.0,
                0,
                3,
                new double[]{0.25, 0.75},
                new int[]{
                        StructuredCoalescentSchedule.SAMPLE,
                        StructuredCoalescentSchedule.COALESCENT
                },
                new int[]{
                        1,
                        StructuredCoalescentSchedule.NO_NODE
                },
                new int[]{
                        StructuredCoalescentSchedule.NO_NODE,
                        0
                },
                new int[]{
                        StructuredCoalescentSchedule.NO_NODE,
                        1
                },
                new int[]{
                        StructuredCoalescentSchedule.NO_NODE,
                        2
                });
    }

    private static int[] filledObservedStates(int nodeCount) {
        int[] states = new int[nodeCount];
        Arrays.fill(states, -1);
        return states;
    }

    private static void assertDoubleArray(double[] expected, double[] actual, double tolerance) {
        assertEquals(expected.length, actual.length);
        for (int i = 0; i < expected.length; i++) {
            assertEquals("entry " + i, expected[i], actual[i], tolerance);
        }
    }
}
