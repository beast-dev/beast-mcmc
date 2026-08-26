/*
 * StructuredCoalescentScheduleTest.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
 */

package test.dr.evomodel.coalescent;

import dr.evolution.io.NewickImporter;
import dr.evomodel.bigfasttree.BigFastTreeIntervals;
import dr.evomodel.bigfasttree.BigFastTreeModel;
import dr.evomodel.coalescent.StructuredCoalescentSchedule;
import dr.evomodel.tree.TreeModel;
import junit.framework.TestCase;

public class StructuredCoalescentScheduleTest extends TestCase {

    public void testScheduleFromTreeIntervals() throws Exception {
        TreeModel tree = fixedTree();
        BigFastTreeIntervals treeIntervals = new BigFastTreeIntervals(tree);

        StructuredCoalescentSchedule schedule = StructuredCoalescentSchedule.fromTreeIntervals(
                tree, treeIntervals, true, false);

        assertEquals(0.0, schedule.initialTime, 0.0);
        assertEquals(1, schedule.initialSampleNode);
        assertEquals(11, schedule.nodeCount);

        assertDoubleArray(new double[]{0.0, 0.5, 0.49, 0.01, 0.5, 0.0, 0.5, 0.5, 0.5, 1.0},
                schedule.intervalLengths, 1.0e-12);
        assertIntArray(new int[]{
                StructuredCoalescentSchedule.SAMPLE,
                StructuredCoalescentSchedule.SAMPLE,
                StructuredCoalescentSchedule.SAMPLE,
                StructuredCoalescentSchedule.COALESCENT,
                StructuredCoalescentSchedule.SAMPLE,
                StructuredCoalescentSchedule.SAMPLE,
                StructuredCoalescentSchedule.COALESCENT,
                StructuredCoalescentSchedule.COALESCENT,
                StructuredCoalescentSchedule.COALESCENT,
                StructuredCoalescentSchedule.COALESCENT
        }, schedule.eventTypes);
        assertIntArray(new int[]{2, 4, 5, -1, 0, 3, -1, -1, -1, -1}, schedule.sampleNodes);
        assertIntArray(new int[]{-1, -1, -1, 1, -1, -1, 0, 4, 7, 8}, schedule.coalescentChild1);
        assertIntArray(new int[]{-1, -1, -1, 2, -1, -1, 6, 5, 3, 9}, schedule.coalescentChild2);
        assertIntArray(new int[]{-1, -1, -1, 6, -1, -1, 7, 9, 8, 10}, schedule.coalescentParent);
    }

    public void testScheduleRejectsBiologicallyInvalidOrder() {
        try {
            new StructuredCoalescentSchedule(
                    0.0,
                    0,
                    3,
                    new double[]{0.1},
                    new int[]{StructuredCoalescentSchedule.COALESCENT},
                    new int[]{StructuredCoalescentSchedule.NO_NODE},
                    new int[]{1},
                    new int[]{2},
                    new int[]{0});
            fail("expected inactive coalescent children to be rejected");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage().contains("inactive children"));
        }
    }

    public void testScheduleRejectsRepeatedSample() {
        try {
            new StructuredCoalescentSchedule(
                    0.0,
                    0,
                    2,
                    new double[]{0.1},
                    new int[]{StructuredCoalescentSchedule.SAMPLE},
                    new int[]{0},
                    new int[]{StructuredCoalescentSchedule.NO_NODE},
                    new int[]{StructuredCoalescentSchedule.NO_NODE},
                    new int[]{StructuredCoalescentSchedule.NO_NODE});
            fail("expected repeated sample to be rejected");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage().contains("already active"));
        }
    }

    public void testScheduleRejectsUnresolvedRoot() {
        try {
            new StructuredCoalescentSchedule(
                    0.0,
                    0,
                    2,
                    new double[]{0.1},
                    new int[]{StructuredCoalescentSchedule.SAMPLE},
                    new int[]{1},
                    new int[]{StructuredCoalescentSchedule.NO_NODE},
                    new int[]{StructuredCoalescentSchedule.NO_NODE},
                    new int[]{StructuredCoalescentSchedule.NO_NODE});
            fail("expected unresolved root to be rejected");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage().contains("ends with 2 active lineages"));
        }
    }

    private static TreeModel fixedTree() throws Exception {
        NewickImporter importer = new NewickImporter(
                "(((0:0.5,(1:1.0,2:1.0)n6:1.0)n7:1.0,3:1.5)n8:1.0," +
                        "(4:2.0,5:1.51)n9:1.5)n10;");
        return new BigFastTreeModel(importer.importTree(null));
    }

    private static void assertDoubleArray(double[] expected, double[] actual, double tolerance) {
        assertEquals(expected.length, actual.length);
        for (int i = 0; i < expected.length; i++) {
            assertEquals("entry " + i, expected[i], actual[i], tolerance);
        }
    }

    private static void assertIntArray(int[] expected, int[] actual) {
        assertEquals(expected.length, actual.length);
        for (int i = 0; i < expected.length; i++) {
            assertEquals("entry " + i, expected[i], actual[i]);
        }
    }
}
