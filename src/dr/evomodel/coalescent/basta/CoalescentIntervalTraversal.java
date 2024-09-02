/*
 * CoalescentIntervalTraversal.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 *
 */

package dr.evomodel.coalescent.basta;

import dr.evolution.coalescent.IntervalType;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.bigfasttree.BigFastTreeIntervals;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.TreeTraversal;

import java.util.*;

import static dr.evomodel.coalescent.basta.ProcessOnCoalescentIntervalDelegate.*;

/**
 * @author Marc A Suchard
 * @author Guy Baele
 */
public class CoalescentIntervalTraversal extends TreeTraversal {

    private final BigFastTreeIntervals treeIntervals;
    private final int numberSubIntervals;

    private int currentMatrixNumber;
    private int currentLikelihoodInterval;

    protected CoalescentIntervalTraversal(final Tree tree,
                                          final BigFastTreeIntervals treeIntervals,
                                          final BranchRateModel branchRateModel,
                                          final int numberSubIntervals) {
        super(tree, branchRateModel, TraversalType.REVERSE_LEVEL_ORDER);

        assert tree instanceof TreeModel;

        this.treeIntervals = treeIntervals;
        this.numberSubIntervals = numberSubIntervals;
    }

    @Override
    public final void dispatchTreeTraversalCollectBranchAndNodeOperations() {
        matrixOperations.clear();

        if (SWAP_API) {
            branchIntervalOperationList.clear();
        } else {
            branchIntervalOperations.clear();
            intervalStarts.clear();
        }

        if (traversalType == TraversalType.REVERSE_LEVEL_ORDER) {
            traverseReverseCoalescentLevelOrder();
        } else {
            assert false : "Unknown traversal type";
        }
    }

    public List<BranchIntervalOperation> getBranchIntervalOperations() {
        if (SWAP_API) {
            return branchIntervalOperationList.getOperations();
        } else {
            return branchIntervalOperations;
        }
    }

    public List<TransitionMatrixOperation> getMatrixOperations() {
        return matrixOperations;
    }

    public int getCoalescentIntervalCount() {
        return currentLikelihoodInterval + 1;
    }

    public List<Integer> getIntervalStarts() {
        if (SWAP_API) {
            return branchIntervalOperationList.getStarts();
        } else {
            return intervalStarts;
        }
    }

    protected final double computeRateScaledIntervalLength(final Tree tree, final NodeRef node, double length) {
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

    static class ActiveNodesForInterval implements Set<NodeRef> {

        private static final boolean DEBUG = true;

        private final Set<NodeRef> activeSet;
        private final int[] currentOffset;
        private final int[] executionOrder;
        private final int stride;
        private final List<NodeRef> intervalNodeOrder;

        public ActiveNodesForInterval(int maximumSize) {
            activeSet = new HashSet<>();
            intervalNodeOrder = new ArrayList<>();
            currentOffset = new int[maximumSize];
            executionOrder = new int[maximumSize];
            stride = maximumSize;
        }

        public Set<NodeRef> copy() {
            return new HashSet<>(activeSet);
        }

        private void test(NodeRef node) {
            if (!activeSet.contains(node)) {
                throw new RuntimeException("Not in active set");
            }
        }

        public int getCurrentOffset(NodeRef node) {
            if (DEBUG) test(node);
            return currentOffset[node.getNumber()];
        }

        public int getActiveBuffer(NodeRef node) {
            if (DEBUG) test(node);
            int currentOffset = getCurrentOffset(node);
            if (currentOffset > 0) {
                ++currentOffset;
            }
            return currentOffset * stride + node.getNumber();
        }

        public int getAccumulationBuffer(NodeRef node) {
            return  stride + node.getNumber();
        }

        public int getExecutionOrder(NodeRef node) {
            if (DEBUG) test(node);
            return executionOrder[node.getNumber()];
        }

        public void incrementActiveBuffer(NodeRef node) {
            if (DEBUG) test(node);
            ++currentOffset[node.getNumber()];
        }

        public void incrementExecutionOrder(NodeRef node) {
            if (DEBUG) test(node);
            ++executionOrder[node.getNumber()];
        }

        public void setExecutionOrder(NodeRef node, int value) {
            if (DEBUG) test(node);
            executionOrder[node.getNumber()] = value;
        }

        public int getNodeOrder(NodeRef node) {
            for (int i = 0; i < intervalNodeOrder.size(); ++i) {
                if (node == intervalNodeOrder.get(i)) return i;
            }
            return -1;
        }

        @Override
        public int size() {
            return activeSet.size();
        }

        @Override
        public boolean isEmpty() {
            return activeSet.isEmpty();
        }

        @Override
        public boolean contains(Object o) {
            return activeSet.contains(o);
        }

        @Override
        public Iterator<NodeRef> iterator() {
            return activeSet.iterator();
        }

        @Override
        public Object[] toArray() {
            throw new UnsupportedOperationException();
        }

        @Override
        public <T> T[] toArray(T[] a) {
            throw new UnsupportedOperationException();
        }

        @Override
        public boolean add(NodeRef node) {
            intervalNodeOrder.add(node);
            return activeSet.add(node);
        }

        @Override
        public boolean remove(Object o) {
            return activeSet.remove(o);
        }

        @Override
        public boolean containsAll(Collection<?> c) {
            throw new UnsupportedOperationException();
        }

        @Override
        public boolean addAll(Collection<? extends NodeRef> c) {
            throw new UnsupportedOperationException();
        }

        @Override
        public boolean retainAll(Collection<?> c) {
            throw new UnsupportedOperationException();
        }

        @Override
        public boolean removeAll(Collection<?> c) {
            throw new UnsupportedOperationException();
        }

        @Override
        public void clear() {
            throw new UnsupportedOperationException();
        }
    }

    //TODO incomplete implementation but also never used
    private int determineStartingInterval() {
        int startingInterval = 0;
        for (int node = 0; node < updateNode.length; ++node) {
            if (updateNode[node]) {

            }
        }
        return startingInterval;
    }

    private void traverseReverseCoalescentLevelOrder() {

        currentLikelihoodInterval = 0;
        currentMatrixNumber = -1;

        // Rebuild active nodes from scratch; TODO cache
        ActiveNodesForInterval activeNodesForInterval = new ActiveNodesForInterval(treeModel.getNodeCount());
        activeNodesForInterval.add(treeIntervals.getSamplingNode(-1)); // Most recent sampled taxon

        intervalStarts.add(0);

        for (int interval = 0; interval < treeIntervals.getIntervalCount(); ++interval) {

            final IntervalType type = treeIntervals.getIntervalType(interval);
            if (type == IntervalType.COALESCENT) {
                processCoalescentEvent(interval, activeNodesForInterval);
            } else if (type == IntervalType.SAMPLE) {
                processSamplingEvent(interval, activeNodesForInterval);
            } else {
                throw new RuntimeException("Unknown interval type");
            }

            if (interval == (treeIntervals.getIntervalCount() - 1)) {
                if (type != IntervalType.COALESCENT) {
                    throw new RuntimeException("Not a coalescence at top");
                }
            }
        }
    }

    @SuppressWarnings("unused")
    private int getDecompositionNumber(NodeRef node) {
        return 0;
    }

    private int computeTransmissionProbabilities(int subInterval, NodeRef node, double length) {

        final int matrixNumber = subInterval; // TODO generalize

        if (matrixNumber != currentMatrixNumber) { // TODO should cache by (decomposition, length)-pair

            final double rateScaledLength = computeRateScaledIntervalLength(treeModel, node, length);

            matrixOperations.add(
                    new TransitionMatrixOperation(
                            matrixNumber,
                            getDecompositionNumber(node),
                            rateScaledLength));

            currentMatrixNumber = matrixNumber;
        }

        return matrixNumber;
    }
    
    private void propagateTransmissionProbabilities(int subInterval, NodeRef node, double length,
                                                    ActiveNodesForInterval activeNodesForInterval) {

        final int inputBuffer1 = activeNodesForInterval.getActiveBuffer(node);
        activeNodesForInterval.incrementActiveBuffer(node);
        final int outputBuffer = activeNodesForInterval.getActiveBuffer(node);
        final int executionOrder = activeNodesForInterval.getExecutionOrder(node) + 1;

        final int inputMatrix1 = computeTransmissionProbabilities(subInterval, node, length);

        BranchIntervalOperation operation = new BranchIntervalOperation(
                outputBuffer,
                inputBuffer1, -1,
                inputMatrix1, -1,
                outputBuffer, -1,
                length, executionOrder, subInterval);

        if (SWAP_API) {
            branchIntervalOperationList.addOperation(operation);
        } else {
            branchIntervalOperations.add(operation);
        }

        activeNodesForInterval.setExecutionOrder(node, executionOrder);
    }

    private void coalescenceTransmissionProbabilities(int subInterval, NodeRef nodeAtTopOfInterval,
                                                      NodeRef leftChild, NodeRef rightChild, double length,
                                                      ActiveNodesForInterval activeNodesForInterval) {
        
        final int inputBuffer1 = activeNodesForInterval.getActiveBuffer(leftChild);
        final int inputBuffer2 = activeNodesForInterval.getActiveBuffer(rightChild);

        final int extraBuffer1 = activeNodesForInterval.getAccumulationBuffer(leftChild);
        final int extraBuffer2 = activeNodesForInterval.getAccumulationBuffer(rightChild);

        final int outputBuffer = activeNodesForInterval.getActiveBuffer(nodeAtTopOfInterval);
        final int executionOrder = Math.max(
                activeNodesForInterval.getExecutionOrder(leftChild),
                activeNodesForInterval.getExecutionOrder(rightChild)) + 1;

        final int inputMatrix1 = computeTransmissionProbabilities(subInterval,leftChild, length);
        final int inputMatrix2 = computeTransmissionProbabilities(subInterval,leftChild, length);

        BranchIntervalOperation operation = new BranchIntervalOperation(
                outputBuffer,
                inputBuffer1, inputBuffer2,
                inputMatrix1, inputMatrix2,
                extraBuffer1, extraBuffer2,
                length, executionOrder, subInterval);

        if (SWAP_API) {
            branchIntervalOperationList.addOperation(operation);
        } else {
            branchIntervalOperations.add(operation);
        }

        activeNodesForInterval.setExecutionOrder(nodeAtTopOfInterval, executionOrder);
    }

    private static final boolean SWAP_API = false;

    private void processCoalescentEvent(int interval, ActiveNodesForInterval activeNodesForInterval) {

        final NodeRef nodeAtTopOfInterval = treeIntervals.getCoalescentNode(interval);
        final double subIntervalLength = treeIntervals.getInterval(interval) / numberSubIntervals;
        final NodeRef leftChild = treeModel.getChild(nodeAtTopOfInterval, 0);
        final NodeRef rightChild = treeModel.getChild(nodeAtTopOfInterval, 1);

        if (subIntervalLength <= 0.0) {
            throw new RuntimeException("Cannot coalescence in <= 0.0 time");
        }

        // Handle initial sub-intervals
        int subInterval = currentLikelihoodInterval * numberSubIntervals;
        for (int i = 0; i < numberSubIntervals - 1; ++i) {
            for (NodeRef activeNode : activeNodesForInterval) {
                propagateTransmissionProbabilities(subInterval, activeNode, subIntervalLength, activeNodesForInterval);
            }
            ++subInterval;
            ++currentLikelihoodInterval;

            if (SWAP_API) {
                branchIntervalOperationList.addStart();
            } else {
                intervalStarts.add(branchIntervalOperations.size());
            }

        }

        // Handle last sub-interval
        activeNodesForInterval.add(nodeAtTopOfInterval);
        coalescenceTransmissionProbabilities(subInterval, nodeAtTopOfInterval, leftChild, rightChild, subIntervalLength,
                activeNodesForInterval);

        //        System.err.println("Added " + nodeAtTopOfInterval.getNumber());

        boolean leftTest = activeNodesForInterval.remove(leftChild);
        boolean rightTest = activeNodesForInterval.remove(rightChild);

        if (!leftTest || !rightTest) {
            throw new RuntimeException("Missing node");
        }

        for (NodeRef activeNode : activeNodesForInterval) {
            if (activeNode != nodeAtTopOfInterval) { // TODO can get rid of check if re-arrange operations above
                propagateTransmissionProbabilities(subInterval, activeNode, subIntervalLength, activeNodesForInterval);
            }
        }

        ++currentLikelihoodInterval;

        if (SWAP_API) {
            branchIntervalOperationList.addStart();
        } else {
            intervalStarts.add(branchIntervalOperations.size());
        }
    }

    private void processSamplingEvent(int interval, ActiveNodesForInterval activeNodesForInterval) {

        final NodeRef nodeAtTopOfInterval = treeIntervals.getSamplingNode(interval);
        final double intervalLength = treeIntervals.getInterval(interval);

        if (intervalLength > 0.0) {

            final double subIntervalLength = intervalLength / numberSubIntervals;

            // Handle all sub-intervals
            int subInterval = currentLikelihoodInterval * numberSubIntervals;
            for (int i = 0; i < numberSubIntervals ; ++i) {
                for (NodeRef activeNode : activeNodesForInterval) {
                    propagateTransmissionProbabilities(subInterval, activeNode, subIntervalLength, activeNodesForInterval);
                }
                ++subInterval;
                ++currentLikelihoodInterval;

                if (SWAP_API) {
                    branchIntervalOperationList.addStart();
                } else {
                    intervalStarts.add(branchIntervalOperations.size());
                }
            }
        }

        activeNodesForInterval.add(nodeAtTopOfInterval);
    }

    private final List<BranchIntervalOperation> branchIntervalOperations = new ArrayList<>();
    private final List<TransitionMatrixOperation> matrixOperations = new ArrayList<>();
    private final List<Integer> intervalStarts = new ArrayList<>();

    private final BranchIntervalOperationList branchIntervalOperationList = null;
}

