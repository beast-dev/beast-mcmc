/*
 * CoalescentIntervalTraversal.java
 *
 * Copyright (c) 2002-2022 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

/**
 * @author Marc A Suchard
 * @author Guy Baele
 */
public class CoalescentIntervalTraversal extends TreeTraversal {

    private final BigFastTreeIntervals treeIntervals;
    private final List<Set<NodeRef>> activeNodesForAllIntervals;
    private final int numberSubIntervals;

    protected CoalescentIntervalTraversal(final Tree tree,
                                          final BigFastTreeIntervals treeIntervals,
                                          final BranchRateModel branchRateModel,
                                          final int numberSubIntervals) {
        super(tree, branchRateModel, TraversalType.REVERSE_LEVEL_ORDER);

        assert tree instanceof TreeModel;

        this.treeIntervals = treeIntervals;
        this.numberSubIntervals = numberSubIntervals;
        this.activeNodesForAllIntervals = new ArrayList<>();
    }

    @Override
    public final void dispatchTreeTraversalCollectBranchAndNodeOperations() {
        branchIntervalOperations.clear();
        otherOperations.clear();

        if (traversalType == TraversalType.REVERSE_LEVEL_ORDER) {
            traverseReverseCoalescentLevelOrder();
        } else {
            assert false : "Unknown traversal type";
        }
    }

    public List<ProcessOnCoalescentIntervalDelegate.BranchIntervalOperation> getBranchIntervalOperations() {
        return branchIntervalOperations;
    }

    public List<ProcessOnCoalescentIntervalDelegate.OtherOperation> getOtherOperations() {
        return otherOperations;
    }

    class ActiveNodesForInterval implements Set<NodeRef> {

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
//            return node.getNumber() * stride + getCurrentOffset(node);
            return 1000 + getCurrentOffset(node) * 1000 + node.getNumber();
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

        int startingInterval = 0;


        // Rebuild active nodes from scratch; TODO cache
        ActiveNodesForInterval activeNodesForInterval = new ActiveNodesForInterval(treeModel.getNodeCount());
        activeNodesForInterval.add(treeIntervals.getSamplingNode(-1)); // Most recent sampled taxon

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

            activeNodesForAllIntervals.add(activeNodesForInterval.copy());
        }

        if (false) {
            for (ProcessOnCoalescentIntervalDelegate.BranchIntervalOperation op : branchIntervalOperations) {
                System.err.println(op);
            }
            System.exit(-1);
        }
    }
    
    private void propagateTransmissionProbabilities(int subInterval, NodeRef node, double length,
                                                    ActiveNodesForInterval activeNodesForInterval) {

        final int inputBuffer1 = activeNodesForInterval.getActiveBuffer(node);
        activeNodesForInterval.incrementActiveBuffer(node);
        final int outputBuffer = activeNodesForInterval.getActiveBuffer(node);
        final int executionOrder = activeNodesForInterval.getExecutionOrder(node) + 1;

        branchIntervalOperations.add(
                new ProcessOnCoalescentIntervalDelegate.BranchIntervalOperation(
                        outputBuffer, inputBuffer1, -1, length, executionOrder, subInterval));

        activeNodesForInterval.setExecutionOrder(node, executionOrder);
    }

    private void coalescenceTransmissionProbabilities(int subInterval, NodeRef nodeAtTopOfInterval,
                                                      NodeRef leftChild, NodeRef rightChild, double length,
                                                      ActiveNodesForInterval activeNodesForInterval) {
        
        final int inputBuffer1 = activeNodesForInterval.getActiveBuffer(leftChild);
        final int inputBuffer2 = activeNodesForInterval.getActiveBuffer(rightChild);

        final int outputBuffer = activeNodesForInterval.getActiveBuffer(nodeAtTopOfInterval);
        final int executionOrder = Math.max(
                activeNodesForInterval.getExecutionOrder(leftChild),
                activeNodesForInterval.getExecutionOrder(rightChild)) + 1;

        branchIntervalOperations.add(
                new ProcessOnCoalescentIntervalDelegate.BranchIntervalOperation(
                        outputBuffer, inputBuffer1, inputBuffer2, length, executionOrder, subInterval));

        activeNodesForInterval.setExecutionOrder(nodeAtTopOfInterval, executionOrder);
    }

    private void processCoalescentEvent(int interval, ActiveNodesForInterval activeNodesForInterval) {

        final NodeRef nodeAtTopOfInterval = treeIntervals.getCoalescentNode(interval);
        final double subIntervalLength = treeIntervals.getInterval(interval) / numberSubIntervals;
        final NodeRef leftChild = treeModel.getChild(nodeAtTopOfInterval, 0);
        final NodeRef rightChild = treeModel.getChild(nodeAtTopOfInterval, 1);

        if (subIntervalLength <= 0.0) {
            throw new RuntimeException("Cannot coalescence in <= 0.0 time");
        }

        // Handle initial sub-intervals
        int subInterval = interval * numberSubIntervals;
        for (int i = 0; i < numberSubIntervals - 1; ++i) {
            for (NodeRef activeNode : activeNodesForInterval) {
                propagateTransmissionProbabilities(subInterval, activeNode, subIntervalLength, activeNodesForInterval);
            }
            ++subInterval;
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
            if (activeNode != nodeAtTopOfInterval) {
                propagateTransmissionProbabilities(subInterval, activeNode, subIntervalLength, activeNodesForInterval);
            }
        }
    }

    private void processSamplingEvent(int interval, ActiveNodesForInterval activeNodesForInterval) {
        final NodeRef nodeAtTopOfInterval = treeIntervals.getSamplingNode(interval);
        final double intervalLength = treeIntervals.getInterval(interval);

        if (intervalLength > 0.0) {

            final double subIntervalLength = intervalLength / numberSubIntervals;

            // Handle all sub-intervals
            int subInterval = interval * numberSubIntervals;
            for (int i = 0; i < numberSubIntervals ; ++i) {
                for (NodeRef activeNode : activeNodesForInterval) {
                    propagateTransmissionProbabilities(subInterval, activeNode, subIntervalLength, activeNodesForInterval);
                }
                ++subInterval;
            }
        }

//        System.err.println("Adding " + nodeAtTopOfInterval.getNumber());
        activeNodesForInterval.add(nodeAtTopOfInterval);
    }

    private final List<ProcessOnCoalescentIntervalDelegate.BranchIntervalOperation> branchIntervalOperations = new ArrayList<>();
    private final List<ProcessOnCoalescentIntervalDelegate.OtherOperation> otherOperations = new ArrayList<>();
}

