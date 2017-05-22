/*
 * SimulationTreeTraversal.java
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

package dr.evomodel.treedatalikelihood;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.BranchRateModel;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by msuchard on 10/6/16.
 */
public class SimulationTreeTraversal extends TreeTraversal {


    public SimulationTreeTraversal(final Tree treeModel,
                                   final BranchRateModel branchRateModel,
                                   final TraversalType traversalType) {
        super(treeModel, branchRateModel, traversalType);
    }

    @Override
    public final void dispatchTreeTraversalCollectBranchAndNodeOperations() {
        branchNodeOperations.clear();
        nodeOperations.clear();

        switch (traversalType) {

            case PRE_ORDER:
                traversePreOrder(treeModel);
                break;
            default:
                assert false : "Unknown traversal type";
        }
    }

    public List<DataLikelihoodDelegate.BranchNodeOperation> getBranchNodeOperations() {
        return branchNodeOperations;
    }

    public List<DataLikelihoodDelegate.NodeOperation> getNodeOperations() {
        return nodeOperations;
    }

    /**
     * Traverse the tree in pre order.
     *
     * @param tree tree
     * @return boolean
     */
    private void traversePreOrder(Tree tree) {
        traversePreOrder(tree, tree.getRoot(), null, null);
    }

    /**
     * Traverse the tree in pre order.
     *
     * @param tree tree
     * @param node node
     * @return boolean
     */
    private boolean traversePreOrder(final Tree tree, final NodeRef node, final NodeRef parent, final NodeRef sibling) {

        boolean update = (parent == null) ? true : false;

        final int nodeNum = node.getNumber();

        // First update the transition probability matrix(ices) for this branch
        if (parent != null && updateNode[nodeNum]) {
            // @todo - at the moment a matrix is updated even if a branch length doesn't change

            branchNodeOperations.add(new DataLikelihoodDelegate.BranchNodeOperation(
                    nodeNum, parent.getNumber(), computeBranchLength(tree, node)
            ));

            nodeOperations.add(new ProcessOnTreeDelegate.NodeOperation(
                    parent.getNumber(),
                    nodeNum,
                    sibling.getNumber()));

            update = true;
        }

        if (!tree.isExternal(node)) {
            final NodeRef child1 = tree.getChild(node, 0);
            final NodeRef child2 = tree.getChild(node, 1);

            if (update) {
                final boolean update1 = traversePreOrder(tree, child1, node, child2);
                final boolean update2 = traversePreOrder(tree, child2, node, child1);
            }
        }

        return update;
    }

    private final List<DataLikelihoodDelegate.BranchNodeOperation> branchNodeOperations =
            new ArrayList<DataLikelihoodDelegate.BranchNodeOperation>();

    private final List<DataLikelihoodDelegate.NodeOperation> nodeOperations =
            new ArrayList<DataLikelihoodDelegate.NodeOperation>();
}
