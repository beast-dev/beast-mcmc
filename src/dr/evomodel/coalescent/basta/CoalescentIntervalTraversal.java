/*
 * TreeTraversal.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.treedatalikelihood.TreeTraversal;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Marc A Suchard
 * @author Guy Baele
 */
public class CoalescentIntervalTraversal extends TreeTraversal {

    protected CoalescentIntervalTraversal(final Tree treeModel,
                                          final BranchRateModel branchRateModel) {
        super(treeModel, branchRateModel, TraversalType.REVERSE_LEVEL_ORDER);
    }

    @Override
    public final void dispatchTreeTraversalCollectBranchAndNodeOperations() {
        branchIntervalOperations.clear();
        otherOperations.clear();

        switch (traversalType) {
            case REVERSE_LEVEL_ORDER:
                traverseReverseCoalescentLevelOrder(treeModel);
                break;
            default:
                assert false : "Unknown traversal type";
        }
    }

    public List<ProcessOnCoalescentIntervalDelegate.BranchIntervalOperation> getBranchIntervalOperations() {
        return branchIntervalOperations;
    }

    public List<ProcessOnCoalescentIntervalDelegate.OtherOperation> getOtherOperations() {
        return otherOperations;
    }


    private void traverseReverseCoalescentLevelOrder(Tree tree) {
        traverseReverseCoalescentLevelOrder(tree, tree.getRoot(), null, null);
    }

    private void traverseReverseCoalescentLevelOrder(Tree tree, NodeRef node1, NodeRef node2, NodeRef node3) {
        // TODO - How does this work?
    }

    List<ProcessOnCoalescentIntervalDelegate.BranchIntervalOperation> branchIntervalOperations = new ArrayList<>();
    List<ProcessOnCoalescentIntervalDelegate.OtherOperation> otherOperations = new ArrayList<>();
}

