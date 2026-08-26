/*
 * CoalescentIntervalTraversalScheduleTest.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
 */

package test.dr.evomodel.coalescent.basta;

import dr.evolution.io.NewickImporter;
import dr.evolution.tree.Tree;
import dr.evomodel.bigfasttree.BigFastTreeIntervals;
import dr.evomodel.bigfasttree.BigFastTreeModel;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
import dr.evomodel.coalescent.basta.CoalescentIntervalTraversal;
import dr.evomodel.coalescent.basta.ProcessOnCoalescentIntervalDelegate.BranchIntervalOperation;
import dr.evomodel.coalescent.basta.ProcessOnCoalescentIntervalDelegate.TransitionMatrixOperation;
import dr.evomodel.tree.TreeModel;
import junit.framework.TestCase;

import java.util.List;

public class CoalescentIntervalTraversalScheduleTest extends TestCase {

    public void testBastaLoweringPreservesCoalescentEventNodes() throws Exception {
        TreeModel tree = fixedTree();
        TestTraversal traversal = new TestTraversal(tree, new BigFastTreeIntervals(tree));

        traversal.dispatchTreeTraversalCollectBranchAndNodeOperations();

        List<BranchIntervalOperation> operations = traversal.getBranchIntervalOperations();
        List<Integer> intervalStarts = traversal.getIntervalStarts();
        List<TransitionMatrixOperation> matrixOperations = traversal.getMatrixOperations();

        assertEquals(Integer.valueOf(0), intervalStarts.get(0));
        assertEquals(Integer.valueOf(operations.size()), intervalStarts.get(intervalStarts.size() - 1));
        for (int i = 1; i < intervalStarts.size(); i++) {
            assertTrue("empty interval group at " + i, intervalStarts.get(i) > intervalStarts.get(i - 1));
        }

        assertDoubleArray(new double[]{0.5, 0.49, 0.01, 0.5, 0.5, 0.5, 0.5, 1.0},
                matrixTimes(matrixOperations), 1.0e-12);
        for (int i = 0; i < matrixOperations.size(); i++) {
            assertEquals(i, matrixOperations.get(i).outputBuffer);
            assertEquals(0, matrixOperations.get(i).decompositionBuffer);
        }

        int[][] expectedCoalescences = {
                {1, 2, 6},
                {0, 6, 7},
                {4, 5, 9},
                {7, 3, 8},
                {8, 9, 10}
        };
        double[] expectedLengths = {0.01, 0.5, 0.5, 0.5, 1.0};

        int coalescentIndex = 0;
        int nodeCount = tree.getNodeCount();
        for (BranchIntervalOperation operation : operations) {
            if (operation.inputBuffer2 < 0) {
                assertEquals(nodeFromBuffer(operation.outputBuffer, nodeCount),
                        nodeFromBuffer(operation.inputBuffer1, nodeCount));
                continue;
            }

            assertTrue("too many coalescent operations", coalescentIndex < expectedCoalescences.length);
            assertEquals(expectedCoalescences[coalescentIndex][0],
                    nodeFromBuffer(operation.inputBuffer1, nodeCount));
            assertEquals(expectedCoalescences[coalescentIndex][1],
                    nodeFromBuffer(operation.inputBuffer2, nodeCount));
            assertEquals(expectedCoalescences[coalescentIndex][2],
                    nodeFromBuffer(operation.outputBuffer, nodeCount));
            assertEquals(expectedLengths[coalescentIndex], operation.intervalLength, 1.0e-12);
            assertEquals(operation.intervalNumber, operation.inputMatrix1);
            assertEquals(operation.intervalNumber, operation.inputMatrix2);
            coalescentIndex++;
        }

        assertEquals(expectedCoalescences.length, coalescentIndex);
    }

    private static double[] matrixTimes(List<TransitionMatrixOperation> matrixOperations) {
        double[] times = new double[matrixOperations.size()];
        for (int i = 0; i < matrixOperations.size(); i++) {
            times[i] = matrixOperations.get(i).time;
        }
        return times;
    }

    private static int nodeFromBuffer(int buffer, int nodeCount) {
        return buffer % nodeCount;
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

    private static final class TestTraversal extends CoalescentIntervalTraversal {
        private TestTraversal(Tree tree, BigFastTreeIntervals treeIntervals) {
            super(tree, treeIntervals, new DefaultBranchRateModel(), 1, true);
        }
    }
}
