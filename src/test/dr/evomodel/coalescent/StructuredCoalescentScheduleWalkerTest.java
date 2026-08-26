/*
 * StructuredCoalescentScheduleWalkerTest.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
 */

package test.dr.evomodel.coalescent;

import dr.evolution.io.NewickImporter;
import dr.evomodel.bigfasttree.BigFastTreeIntervals;
import dr.evomodel.bigfasttree.BigFastTreeModel;
import dr.evomodel.coalescent.StructuredCoalescentActiveLineages;
import dr.evomodel.coalescent.StructuredCoalescentSchedule;
import dr.evomodel.coalescent.StructuredCoalescentScheduleWalker;
import dr.evomodel.tree.TreeModel;
import junit.framework.TestCase;

import java.util.ArrayList;
import java.util.List;

public class StructuredCoalescentScheduleWalkerTest extends TestCase {

    public void testWalkerOwnsBiologicalActiveLineageOrder() throws Exception {
        TreeModel tree = fixedTree();
        StructuredCoalescentSchedule schedule = StructuredCoalescentSchedule.fromTreeIntervals(
                tree, new BigFastTreeIntervals(tree), true, false);
        RecordingVisitor visitor = new RecordingVisitor();

        StructuredCoalescentScheduleWalker.walk(schedule, visitor);

        assertStringArray(new String[]{
                "start 1 [1]",
                "interval 0 1 [1]",
                "sample 0 2 [1]",
                "after 0 1 [1,2]",
                "interval 1 1 [1,2]",
                "sample 1 4 [1,2]",
                "after 1 1 [1,2,4]",
                "interval 2 1 [1,2,4]",
                "sample 2 5 [1,2,4]",
                "after 2 1 [1,2,4,5]",
                "interval 3 2 [1,2,4,5]",
                "coalescent 3 1 2 6 [1,2,4,5]",
                "after 3 2 [4,5,6]",
                "interval 4 1 [4,5,6]",
                "sample 4 0 [4,5,6]",
                "after 4 1 [4,5,6,0]",
                "interval 5 1 [4,5,6,0]",
                "sample 5 3 [4,5,6,0]",
                "after 5 1 [4,5,6,0,3]",
                "interval 6 2 [4,5,6,0,3]",
                "coalescent 6 0 6 7 [4,5,6,0,3]",
                "after 6 2 [4,5,3,7]",
                "interval 7 2 [4,5,3,7]",
                "coalescent 7 4 5 9 [4,5,3,7]",
                "after 7 2 [3,7,9]",
                "interval 8 2 [3,7,9]",
                "coalescent 8 7 3 8 [3,7,9]",
                "after 8 2 [9,8]",
                "interval 9 2 [9,8]",
                "coalescent 9 8 9 10 [9,8]",
                "after 9 2 [10]",
                "finish [10]"
        }, visitor.records);
    }

    private static TreeModel fixedTree() throws Exception {
        NewickImporter importer = new NewickImporter(
                "(((0:0.5,(1:1.0,2:1.0)n6:1.0)n7:1.0,3:1.5)n8:1.0," +
                        "(4:2.0,5:1.51)n9:1.5)n10;");
        return new BigFastTreeModel(importer.importTree(null));
    }

    private static void assertStringArray(String[] expected, List<String> actual) {
        assertEquals(expected.length, actual.size());
        for (int i = 0; i < expected.length; i++) {
            assertEquals("entry " + i, expected[i], actual.get(i));
        }
    }

    private static String active(StructuredCoalescentActiveLineages activeLineages) {
        StringBuilder builder = new StringBuilder("[");
        for (int i = 0; i < activeLineages.getActiveCount(); i++) {
            if (i > 0) {
                builder.append(',');
            }
            builder.append(activeLineages.getActiveNode(i));
        }
        builder.append(']');
        return builder.toString();
    }

    private static final class RecordingVisitor extends StructuredCoalescentScheduleWalker.Adapter {
        private final List<String> records = new ArrayList<>();

        @Override
        public void initialSample(int sampleNode, StructuredCoalescentActiveLineages activeLineages) {
            records.add("start " + sampleNode + " " + active(activeLineages));
        }

        @Override
        public void interval(int interval, double intervalLength, int eventType,
                             StructuredCoalescentActiveLineages activeLineages) {
            records.add("interval " + interval + " " + eventType + " " + active(activeLineages));
        }

        @Override
        public void sampleEvent(int interval, double intervalLength, int sampleNode,
                                StructuredCoalescentActiveLineages activeLineages) {
            records.add("sample " + interval + " " + sampleNode + " " + active(activeLineages));
        }

        @Override
        public void coalescentEvent(int interval, double intervalLength, int leftChild, int rightChild, int parent,
                                    StructuredCoalescentActiveLineages activeLineages) {
            records.add("coalescent " + interval + " " + leftChild + " " + rightChild + " " + parent + " " +
                    active(activeLineages));
        }

        @Override
        public void afterEvent(int interval, int eventType, StructuredCoalescentActiveLineages activeLineages) {
            records.add("after " + interval + " " + eventType + " " + active(activeLineages));
        }

        @Override
        public void finish(StructuredCoalescentActiveLineages activeLineages) {
            records.add("finish " + active(activeLineages));
        }
    }
}
