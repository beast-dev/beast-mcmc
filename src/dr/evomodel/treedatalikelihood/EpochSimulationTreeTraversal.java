/*
 * EpochSimulationTreeTraversal.java
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

package dr.evomodel.treedatalikelihood;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.BranchRateModel;

/**
 * Pre-order traversal for an epoch model, over the tree that is augmented with a degree-2 node at each epoch
 * transition time that a branch spans (see EpochLikelihoodTraversal, which assigns the degree-2 nodes when the
 * likelihood is evaluated; this traversal only reads them).
 * <p>
 * For the branch above node n, with augmented nodes a_0, ..., a_(k-1) from the bottom to the top, and parent p, the
 * operations go from the top to the bottom:
 * <ul>
 * <li>NodeOperation (p; a_(k-1); sibling of n) and then (a_(j+1); a_j; NO_CHILD) and (a_0; n; NO_CHILD), and</li>
 * <li>a BranchNodeOperation for each segment.</li>
 * </ul>
 * The sibling is read as its top augmented node. Each NodeOperation records the level of the node that it
 * computes, its depth below the root of the augmented tree.
 * <p>
 * Like SimulationTreeTraversal, an operation is added for each node that is marked for updating, and the children of
 * a node are only visited when the node is marked. The rate of each branch is assumed constant over its segments.
 *
 * @author Marc A Suchard
 */
public class EpochSimulationTreeTraversal extends SimulationTreeTraversal {

    private final AugmentedNodeRegistry registry;
    private final EpochEvolutionaryProcessDelegate epochDelegate;

    private final int[] chain;

    public EpochSimulationTreeTraversal(final Tree treeModel,
                                        final BranchRateModel branchRateModel,
                                        final TraversalType traversalType,
                                        final EpochEvolutionaryProcessDelegate epochDelegate) {
        super(treeModel, branchRateModel, traversalType);

        if (traversalType != TraversalType.PRE_ORDER) {
            throw new IllegalArgumentException("Unknown traversal type");
        }

        this.epochDelegate = epochDelegate;
        this.registry = epochDelegate.getAugmentedNodeRegistry();
        this.chain = new int[registry.getBoundaryCount()];
    }

    @Override
    public void dispatchTreeTraversalCollectBranchAndNodeOperations() {
        branchNodeOperations.clear();
        nodeOperations.clear();

        if (!registry.hasBoundaries(epochDelegate.getEpochTransitionTimes())) {
            throw new IllegalStateException("The degree-2 nodes are out of date; evaluate the likelihood first");
        }

        traversePreOrder(treeModel, treeModel.getRoot(), 0);
    }

    /**
     * @param level depth of node below the root
     */
    private void traversePreOrder(final Tree tree, final NodeRef node, final int level) {
        if (tree.isExternal(node)) {
            return;
        }

        final NodeRef child1 = tree.getChild(node, 0);
        final NodeRef child2 = tree.getChild(node, 1);

        int childLevel = addOperations(tree, child1, node, child2, level);
        if (childLevel >= 0) {
            traversePreOrder(tree, child1, childLevel);
        }

        childLevel = addOperations(tree, child2, node, child1, level);
        if (childLevel >= 0) {
            traversePreOrder(tree, child2, childLevel);
        }
    }

    /**
     * Adds the operations for the branch above node, from the top to the bottom.
     *
     * @return the level of node, or -1 if the branch is not marked for updating
     */
    private int addOperations(final Tree tree, final NodeRef node, final NodeRef parent, final NodeRef sibling,
                              final int parentLevel) {

        final int nodeNum = node.getNumber();
        if (!updateNode[nodeNum]) {
            return -1;
        }

        final double branchRate;
        synchronized (branchRateModel) {
            branchRate = branchRateModel.getBranchRate(tree, node);
        }

        final int count = registry.copyChain(nodeNum, chain);
        final int first = registry.getFirstBoundary(nodeNum);

        // above the first (highest) segment the sibling meets the branch; the rest have one child
        int siblingInput = registry.getChainTop(sibling.getNumber());

        int upperNode = parent.getNumber();
        double upperHeight = tree.getNodeHeight(parent);
        int level = parentLevel;

        for (int j = count - 1; j >= -1; --j) {
            final int lowerNode = (j >= 0) ? chain[j] : nodeNum;
            final double lowerHeight = (j >= 0) ? registry.getBoundary(first + j) : tree.getNodeHeight(node);
            ++level;

            final double segmentLength = branchRate * (upperHeight - lowerHeight);
            assert segmentLength >= 0.0 : "Negative segment length: " + segmentLength + " for node " + nodeNum;

            branchNodeOperations.add(new DataLikelihoodDelegate.BranchNodeOperation(lowerNode, upperNode,
                    segmentLength));
            nodeOperations.add(new ProcessOnTreeDelegate.NodeOperation(upperNode, lowerNode, siblingInput, level));

            siblingInput = ProcessOnTreeDelegate.NodeOperation.NO_CHILD;
            upperNode = lowerNode;
            upperHeight = lowerHeight;
        }

        return level;
    }
}
