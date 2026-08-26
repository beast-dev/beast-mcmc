/*
 * BastaBranchIntervalOperationBuilder.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
 */

package dr.evomodel.coalescent.basta;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.coalescent.StructuredCoalescentActiveLineages;
import dr.evomodel.coalescent.StructuredCoalescentSchedule;
import dr.evomodel.coalescent.StructuredCoalescentScheduleWalker;

import java.util.Arrays;
import java.util.List;

import static dr.evomodel.coalescent.basta.ProcessOnCoalescentIntervalDelegate.BranchIntervalOperation;
import static dr.evomodel.coalescent.basta.ProcessOnCoalescentIntervalDelegate.TransitionMatrixOperation;

/**
 * Lowers a shared biological structured-coalescent schedule into the existing
 * BASTA operation stream.
 */
final class BastaBranchIntervalOperationBuilder extends StructuredCoalescentScheduleWalker.Adapter {

    private final Tree tree;
    private final BranchRateModel branchRateModel;
    private final int numberSubIntervals;
    private final boolean checkForZeroLengthIntervals;
    private final List<BranchIntervalOperation> branchIntervalOperations;
    private final List<TransitionMatrixOperation> matrixOperations;
    private final List<Integer> intervalStarts;
    private final ActiveBuffers activeBuffers;

    private int currentMatrixNumber;
    private int currentLikelihoodInterval;

    BastaBranchIntervalOperationBuilder(Tree tree,
                                        BranchRateModel branchRateModel,
                                        int numberSubIntervals,
                                        boolean checkForZeroLengthIntervals,
                                        List<BranchIntervalOperation> branchIntervalOperations,
                                        List<TransitionMatrixOperation> matrixOperations,
                                        List<Integer> intervalStarts) {
        if (numberSubIntervals <= 0) {
            throw new IllegalArgumentException("numberSubIntervals must be positive");
        }
        this.tree = tree;
        this.branchRateModel = branchRateModel;
        this.numberSubIntervals = numberSubIntervals;
        this.checkForZeroLengthIntervals = checkForZeroLengthIntervals;
        this.branchIntervalOperations = branchIntervalOperations;
        this.matrixOperations = matrixOperations;
        this.intervalStarts = intervalStarts;
        this.activeBuffers = new ActiveBuffers(tree.getNodeCount());
    }

    void build(StructuredCoalescentSchedule schedule) {
        currentLikelihoodInterval = 0;
        currentMatrixNumber = -1;
        activeBuffers.clear();
        intervalStarts.add(0);
        StructuredCoalescentScheduleWalker.walk(schedule, this);
    }

    int getCoalescentIntervalCount() {
        return currentLikelihoodInterval + 1;
    }

    @Override
    public void sampleEvent(int interval, double intervalLength, int sampleNode,
                            StructuredCoalescentActiveLineages activeLineages) {
        if (intervalLength > 0.0) {
            int subInterval = currentLikelihoodInterval * numberSubIntervals;
            double subIntervalLength = intervalLength / numberSubIntervals;
            for (int i = 0; i < numberSubIntervals; i++) {
                propagateActiveLineages(subInterval, subIntervalLength, activeLineages,
                        StructuredCoalescentSchedule.NO_NODE, StructuredCoalescentSchedule.NO_NODE);
                subInterval++;
                completeLikelihoodInterval();
            }
        }
    }

    @Override
    public void coalescentEvent(int interval, double intervalLength, int leftChild, int rightChild, int parent,
                                StructuredCoalescentActiveLineages activeLineages) {
        double subIntervalLength = intervalLength / numberSubIntervals;
        if (checkForZeroLengthIntervals && subIntervalLength <= 0.0) {
            throw new RuntimeException("Cannot coalesce in <= 0.0 time");
        }

        int subInterval = currentLikelihoodInterval * numberSubIntervals;
        for (int i = 0; i < numberSubIntervals - 1; i++) {
            propagateActiveLineages(subInterval, subIntervalLength, activeLineages,
                    StructuredCoalescentSchedule.NO_NODE, StructuredCoalescentSchedule.NO_NODE);
            subInterval++;
            completeLikelihoodInterval();
        }

        coalescenceTransmissionProbabilities(subInterval, parent, leftChild, rightChild, subIntervalLength);
        propagateActiveLineages(subInterval, subIntervalLength, activeLineages, leftChild, rightChild);
        completeLikelihoodInterval();
    }

    private void propagateActiveLineages(int subInterval, double length,
                                         StructuredCoalescentActiveLineages activeLineages,
                                         int skipNode1, int skipNode2) {
        for (int i = 0; i < activeLineages.getActiveCount(); i++) {
            int activeNode = activeLineages.getActiveNode(i);
            if (activeNode != skipNode1 && activeNode != skipNode2) {
                propagateTransmissionProbabilities(subInterval, activeNode, length);
            }
        }
    }

    private void completeLikelihoodInterval() {
        currentLikelihoodInterval++;
        intervalStarts.add(branchIntervalOperations.size());
    }

    private int getDecompositionNumber(int node) {
        return 0;
    }

    private int computeTransmissionProbabilities(int subInterval, int node, double length) {
        final int matrixNumber = subInterval; // TODO generalize

        if (matrixNumber != currentMatrixNumber) { // TODO should cache by (decomposition, length)-pair
            final double rateScaledLength = computeRateScaledIntervalLength(node, length);

            matrixOperations.add(
                    new TransitionMatrixOperation(
                            matrixNumber,
                            getDecompositionNumber(node),
                            rateScaledLength));

            currentMatrixNumber = matrixNumber;
        }

        return matrixNumber;
    }

    private double computeRateScaledIntervalLength(int nodeNumber, double length) {
        NodeRef node = tree.getNode(nodeNumber);
        final double branchRate;

        synchronized (branchRateModel) {
            branchRate = branchRateModel.getBranchRate(tree, node);
        }

        final double branchIntervalLength = branchRate * length;

        assert branchIntervalLength >= 0.0 : "Negative interval length: " + branchIntervalLength + " for node " +
                node.getNumber() + (tree.isExternal(node) ?
                " (" + tree.getNodeTaxon(node).getId() + ")" : "");

        return branchIntervalLength;
    }

    private void propagateTransmissionProbabilities(int subInterval, int node, double length) {
        final int inputBuffer1 = activeBuffers.getActiveBuffer(node);
        activeBuffers.incrementActiveBuffer(node);
        final int outputBuffer = activeBuffers.getActiveBuffer(node);
        final int executionOrder = activeBuffers.getExecutionOrder(node) + 1;

        final int inputMatrix1 = computeTransmissionProbabilities(subInterval, node, length);

        BranchIntervalOperation operation = new BranchIntervalOperation(
                outputBuffer,
                inputBuffer1, -1,
                inputMatrix1, -1,
                outputBuffer, -1,
                length, executionOrder, subInterval);

        branchIntervalOperations.add(operation);
        activeBuffers.setExecutionOrder(node, executionOrder);
    }

    private void coalescenceTransmissionProbabilities(int subInterval, int parent,
                                                      int leftChild, int rightChild, double length) {
        final int inputBuffer1 = activeBuffers.getActiveBuffer(leftChild);
        final int inputBuffer2 = activeBuffers.getActiveBuffer(rightChild);

        final int extraBuffer1 = activeBuffers.getAccumulationBuffer(leftChild);
        final int extraBuffer2 = activeBuffers.getAccumulationBuffer(rightChild);

        final int outputBuffer = activeBuffers.getActiveBuffer(parent);
        final int executionOrder = Math.max(
                activeBuffers.getExecutionOrder(leftChild),
                activeBuffers.getExecutionOrder(rightChild)) + 1;

        final int inputMatrix1 = computeTransmissionProbabilities(subInterval, leftChild, length);
        final int inputMatrix2 = computeTransmissionProbabilities(subInterval, rightChild, length);

        BranchIntervalOperation operation = new BranchIntervalOperation(
                outputBuffer,
                inputBuffer1, inputBuffer2,
                inputMatrix1, inputMatrix2,
                extraBuffer1, extraBuffer2,
                length, executionOrder, subInterval);

        branchIntervalOperations.add(operation);
        activeBuffers.setExecutionOrder(parent, executionOrder);
    }

    private static final class ActiveBuffers {

        private final int[] currentOffset;
        private final int[] executionOrder;
        private final int stride;

        private ActiveBuffers(int maximumSize) {
            currentOffset = new int[maximumSize];
            executionOrder = new int[maximumSize];
            stride = maximumSize;
        }

        private void clear() {
            Arrays.fill(currentOffset, 0);
            Arrays.fill(executionOrder, 0);
        }

        private int getActiveBuffer(int node) {
            int offset = currentOffset[node];
            if (offset > 0) {
                offset++;
            }
            return offset * stride + node;
        }

        private int getAccumulationBuffer(int node) {
            return stride + node;
        }

        private int getExecutionOrder(int node) {
            return executionOrder[node];
        }

        private void incrementActiveBuffer(int node) {
            currentOffset[node]++;
        }

        private void setExecutionOrder(int node, int value) {
            executionOrder[node] = value;
        }
    }
}
