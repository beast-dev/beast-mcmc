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
import dr.evomodel.coalescent.TreeIntervals;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.TreeTraversal;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * @author Marc A Suchard
 * @author Guy Baele
 */
public class CoalescentIntervalTraversal extends TreeTraversal {

    private final BigFastTreeIntervals treeIntervals;
    private final List<Set<NodeRef>> activeNodesForAllIntervals;

    protected CoalescentIntervalTraversal(final Tree tree,
                                          final BranchRateModel branchRateModel) {
        super(tree, branchRateModel, TraversalType.REVERSE_LEVEL_ORDER);

        assert tree instanceof TreeModel;
        treeIntervals = new BigFastTreeIntervals((TreeModel) tree);

        activeNodesForAllIntervals = new ArrayList<>();
        
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


    
    private void traverseReverseCoalescentLevelOrder() {

        // Rebuild active nodes from scratch; TODO cache
//        activeNodesForAllIntervals.clear();
//        Set<NodeRef> activeNodesForPreviousInterval = new HashSet<>();
        Set<NodeRef> activeNodesForInterval = new HashSet<>();

        for (int interval = 0; interval < treeIntervals.getIntervalCount(); ++interval) {

//            Set<NodeRef> activeNodesForThisInterval =

            final IntervalType type = treeIntervals.getIntervalType(interval);
            if (type == IntervalType.COALESCENT) {
                processCoalescentEvent(interval, activeNodesForInterval);
            } else if (type == IntervalType.SAMPLE) {
                processSamplingEvent(interval, activeNodesForInterval);
            } else {
                throw new RuntimeException("Unknown interval type");
            }
        }
    }

    private void processCoalescentEvent(int interval, Set<NodeRef> activeNodesForInterval) {
        final NodeRef node = treeIntervals.getCoalescentNode(interval);
        final NodeRef leftChild = treeModel.getChild(node, 0);
        final NodeRef rightChild = treeModel.getChild(node, 1);

        boolean leftTest = activeNodesForInterval.remove(leftChild);
        boolean rightTest = activeNodesForInterval.remove(rightChild);

        if (!leftTest || !rightTest) {
            throw new RuntimeException("Missing node");
        }

        activeNodesForInterval.add(node);
    }

    private void processSamplingEvent(int interval, Set<NodeRef> activeNodesForInterval) {
//        NodeRef node = treeIntervals.getSamplingNode(interval);
        final NodeRef node = null;

        activeNodesForInterval.add(node);
    }

    private final List<ProcessOnCoalescentIntervalDelegate.BranchIntervalOperation> branchIntervalOperations = new ArrayList<>();
    private final List<ProcessOnCoalescentIntervalDelegate.OtherOperation> otherOperations = new ArrayList<>();


}

