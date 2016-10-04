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

package dr.evomodel.treedatalikelihood;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.BranchRateModel;

import java.util.*;

/**
 * @author Andrew Rambaut
 * @author Marc A Suchard
 */
public class TreeTraversal {

    public TreeTraversal(final Tree treeModel, final boolean[] updateNode,
                         final BranchRateModel branchRateModel,
                         final TreeDataLikelihood.TraversalType traversalType) {

        this.treeModel = treeModel;
        this.branchRateModel = branchRateModel;
        this.updateNode = updateNode;
        this.traversalType = traversalType;
    }


    public void dispatchTreeTraversalCollectBranchAndNodeOperations() {
        branchOperations.clear();
        nodeOperations.clear();

        switch (traversalType) {

            case POST_ORDER:
                traversePostOrder(treeModel);
                break;
            case REVERSE_LEVEL_ORDER:
                traverseReverseLevelOrder(treeModel);
                break;
            default:
                assert false : "Unknown traversal type";
        }
    }

    /**
     * Traverse the tree in post order.
     *
     * @param tree tree
     * @return boolean
     */
    private void traversePostOrder(Tree tree) {
        traversePostOrder(tree, tree.getRoot());
    }

    /**
     * Traverse the tree in post order.
     *
     * @param tree tree
     * @param node node
     * @return boolean
     */
    private boolean traversePostOrder(Tree tree, NodeRef node) {

        boolean update = false;

        int nodeNum = node.getNumber();

        // First update the transition probability matrix(ices) for this branch
        if (tree.getParent(node) != null && updateNode[nodeNum]) {
            // @todo - at the moment a matrix is updated even if a branch length doesn't change

            addBranchUpdateOperation(tree, node);

            update = true;
        }

        // If the node is internal, update the partial likelihoods.
        if (!tree.isExternal(node)) {

            // Traverse down the two child nodes
            NodeRef child1 = tree.getChild(node, 0);
            final boolean update1 = traversePostOrder(tree, child1);

            NodeRef child2 = tree.getChild(node, 1);
            final boolean update2 = traversePostOrder(tree, child2);

            // If either child node was updated then update this node too
            if (update1 || update2) {

                nodeOperations.add(new DataLikelihoodDelegate.NodeOperation(nodeNum, child1.getNumber(), child2.getNumber()));

                update = true;

            }
        }

        return update;

    }

    /**
     * Traverse the tree in reverse level order.
     *
     * @param tree tree
     */
    private void traverseReverseLevelOrder(final Tree tree) {

        // create a map of all the operations at each particular level
        Map<Integer, List<DataLikelihoodDelegate.NodeOperation>> operationMap =
                new HashMap<Integer, List<DataLikelihoodDelegate.NodeOperation>>();

        traverseLevelOrder(tree, tree.getRoot(), 0, operationMap);

        // get the levels as keys in reverse order (they are currently largest towards
        // the tips) and add the operations to the nodeOperation array.
        List<Integer> keyList = new ArrayList<Integer>(operationMap.keySet());
        Collections.sort(keyList, Collections.reverseOrder());

        for (Integer key : keyList) {
            List<DataLikelihoodDelegate.NodeOperation> opList = operationMap.get(key);
            for (DataLikelihoodDelegate.NodeOperation op : opList) {
                nodeOperations.add(op);
            }
        }
    }

    /**
     * Traverse the tree in level order.
     *
     * @param tree tree
     * @param node node
     * @return boolean
     */
    private boolean traverseLevelOrder(final Tree tree, final NodeRef node,
                                       final int level,
                                       Map<Integer, List<DataLikelihoodDelegate.NodeOperation>> operationMap) {
        boolean update = false;

        int nodeNum = node.getNumber();

        // First update the transition probability matrix(ices) for this branch
        if (tree.getParent(node) != null && updateNode[nodeNum]) {
            // @todo - at the moment a matrix is updated even if a branch length doesn't change

            addBranchUpdateOperation(tree, node);

            update = true;
        }

        // If the node is internal, update the partial likelihoods.
        if (!tree.isExternal(node)) {

            // Traverse down the two child nodes incrementing the level (this will give
            // level order but we will reverse these later
            NodeRef child1 = tree.getChild(node, 0);
            final boolean update1 = traverseLevelOrder(tree, child1, level + 1, operationMap);

            NodeRef child2 = tree.getChild(node, 1);
            final boolean update2 = traverseLevelOrder(tree, child2, level + 1, operationMap);

            // If either child node was updated then update this node too
            if (update1 || update2) {

                List<DataLikelihoodDelegate.NodeOperation> ops = operationMap.get(level);
                if (ops == null) {
                    ops = new ArrayList<DataLikelihoodDelegate.NodeOperation>();
                    operationMap.put(level, ops);
                }
                ops.add(new DataLikelihoodDelegate.NodeOperation(nodeNum, child1.getNumber(), child2.getNumber()));

                update = true;

            }
        }

        return update;
    }

    /**
     * Add this node to the branchOperation list for updating of the transition probability matrix.
     *
     * @param tree tree
     * @param node node
     */
    private void addBranchUpdateOperation(Tree tree, NodeRef node) {
        final double branchRate;

        synchronized (branchRateModel) {
            branchRate = branchRateModel.getBranchRate(tree, node);
        }
        final double parentHeight = tree.getNodeHeight(tree.getParent(node));
        final double nodeHeight = tree.getNodeHeight(node);

        // Get the operational time of the branch
        final double branchLength = branchRate * (parentHeight - nodeHeight);

        assert branchLength > 0.0 : "Negative branch length: " + branchLength + " for node " +
                node.getNumber() + (tree.isExternal(node) ?
                " (" + tree.getNodeTaxon(node).getId() + ")" : "");

        branchOperations.add(new DataLikelihoodDelegate.BranchOperation(node.getNumber(), branchLength));
    }

    public List<DataLikelihoodDelegate.BranchOperation> getBranchOperations() {
        return branchOperations;
    }

    public List<DataLikelihoodDelegate.NodeOperation> getNodeOperations() {
        return nodeOperations;
    }

    final private Tree treeModel;
    final private BranchRateModel branchRateModel;
    final private boolean[] updateNode;

    private final TreeDataLikelihood.TraversalType traversalType;


    private List<DataLikelihoodDelegate.BranchOperation> branchOperations = new ArrayList<DataLikelihoodDelegate.BranchOperation>();
    private List<DataLikelihoodDelegate.NodeOperation> nodeOperations = new ArrayList<DataLikelihoodDelegate.NodeOperation>();

}
