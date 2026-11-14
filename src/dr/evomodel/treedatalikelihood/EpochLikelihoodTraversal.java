/*
 * EpochLikelihoodTraversal.java
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

import java.util.Arrays;

/**
 * Post-order traversal for an epoch model. In addition to the branch and node operations of
 * LikelihoodTreeTraversal, a degree-2 (single-child) node is added to a branch at each epoch transition time that
 * the branch spans, so that each segment of the branch lies within a single epoch.
 * <p>
 * For the branch above node n that crosses transition times t_f, ..., t_(f+k-1), with a_j the augmented node at t_(f+j),
 * <ul>
 * <li>one BranchOperation for each of n, a_0, ..., a_(k-1) updates the matrix of the segment above it,</li>
 * <li>a NodeOperation (a_j; input a_(j-1) or n; NO_CHILD) computes the partials of a_j, and</li>
 * <li>the parent of n reads a_(k-1) instead of n.</li>
 * </ul>
 * Each NodeOperation records its level, the depth below the root of the augmented tree. Operations at the same
 * level are independent. For REVERSE_LEVEL_ORDER the node operations are sorted by decreasing level.
 * The rate of each branch is assumed constant over its segments.
 *
 * @author Marc A Suchard
 */
public class EpochLikelihoodTraversal extends LikelihoodTreeTraversal {

    private final AugmentedNodeRegistry registry;
    private final EpochEvolutionaryProcessDelegate epochDelegate;

    private final int[] chain;
    private double[] boundaries;

    // branches whose chains grow, assigned after every other branch has released (see reassignChains)
    private final int[] pendingNode;
    private final int[] pendingFirst;
    private final int[] pendingCount;

    public EpochLikelihoodTraversal(final Tree treeModel,
                                    final BranchRateModel branchRateModel,
                                    final TraversalType traversalType,
                                    final EpochEvolutionaryProcessDelegate epochDelegate) {
        super(treeModel, branchRateModel, traversalType);

        if (traversalType != TraversalType.POST_ORDER && traversalType != TraversalType.REVERSE_LEVEL_ORDER) {
            throw new IllegalArgumentException("Unknown traversal type");
        }

        this.epochDelegate = epochDelegate;
        this.registry = epochDelegate.getAugmentedNodeRegistry();
        this.chain = new int[registry.getBoundaryCount()];

        final int nodeCount = treeModel.getNodeCount();
        this.pendingNode = new int[nodeCount];
        this.pendingFirst = new int[nodeCount];
        this.pendingCount = new int[nodeCount];
    }

    @Override
    public void dispatchTreeTraversalCollectBranchAndNodeOperations() {
        branchOperations.clear();
        nodeOperations.clear();

        boundaries = epochDelegate.getEpochTransitionTimes();
        if (registry.setBoundaries(boundaries)) {
            updateAllNodes();
        }

        reassignChains(treeModel);

        traversePostOrder(treeModel, treeModel.getRoot(), 0);

        if (traversalType == TraversalType.REVERSE_LEVEL_ORDER) {
            sortNodeOperationsByReverseLevel();
        }
    }

    @Override
    protected int getChildLevel(final Tree tree, final NodeRef child, final int parentLevel) {
        return parentLevel + 1 + registry.getChainLength(child.getNumber());
    }

    @Override
    protected void addBranchOperations(final Tree tree, final NodeRef node) {
        final int nodeNum = node.getNumber();

        final double branchRate;
        synchronized (branchRateModel) {
            branchRate = branchRateModel.getBranchRate(tree, node);
        }

        final int count = registry.copyChain(nodeNum, chain);
        final int first = registry.getFirstBoundary(nodeNum);
        final double upperHeight = tree.getNodeHeight(tree.getParent(node));

        int lowerNode = nodeNum;
        double lowerHeight = tree.getNodeHeight(node);

        for (int j = 0; j <= count; ++j) {
            final double segmentTop = (j < count) ? boundaries[first + j] : upperHeight;

            final double segmentLength = branchRate * (segmentTop - lowerHeight);
            assert segmentLength >= 0.0 : "Negative segment length: " + segmentLength + " for node " + nodeNum;

            branchOperations.add(new DataLikelihoodDelegate.BranchOperation(lowerNode, segmentLength));

            if (j < count) {
                lowerNode = chain[j];
            }
            lowerHeight = segmentTop;
        }
    }

    @Override
    protected int getInputNodeNumber(final NodeRef child) {
        return registry.getChainTop(child.getNumber());
    }

    @Override
    protected void afterNodeOperation(final Tree tree, final NodeRef node, final int level) {
        final int nodeNum = node.getNumber();
        final int count = registry.copyChain(nodeNum, chain);

        int lowerNode = nodeNum;
        for (int j = 0; j < count; ++j) {
            nodeOperations.add(new DataLikelihoodDelegate.NodeOperation(chain[j], lowerNode,
                    DataLikelihoodDelegate.NodeOperation.NO_CHILD, level - 1 - j));
            lowerNode = chain[j];
        }
    }

    /**
     * Gives the branch above each updated node a chain of augmented nodes, one for each transition time strictly
     * between the height of the node and the height of its parent. Chains that shrink or keep their length are
     * assigned first and chains that grow afterwards, so the number in use never exceeds the larger of the old and
     * the new total, and taxonCount * boundaryCount always suffices. The root has no branch, so a node that became
     * the root releases its chain.
     */
    private void reassignChains(final Tree tree) {
        final int root = tree.getRoot().getNumber();
        if (registry.getChainLength(root) > 0) {
            registry.assign(root, 0, 0);
        }

        int growing = 0;
        for (int n = 0; n < tree.getNodeCount(); ++n) {
            if (n == root || !updateNode[n]) {
                continue;
            }
            final NodeRef node = tree.getNode(n);
            assert node.getNumber() == n;

            final int first = firstIndexAbove(boundaries, tree.getNodeHeight(node));
            final int count = Math.max(0,
                    firstIndexAtOrAbove(boundaries, tree.getNodeHeight(tree.getParent(node))) - first);

            if (count <= registry.getChainLength(n)) {
                registry.assign(n, first, count); // releases or relabels only
            } else {
                pendingNode[growing] = n;
                pendingFirst[growing] = first;
                pendingCount[growing] = count;
                ++growing;
            }
        }

        for (int k = 0; k < growing; ++k) { // takes only
            registry.assign(pendingNode[k], pendingFirst[k], pendingCount[k]);
        }
    }

    /**
     * @return the smallest index i with times[i] > value, or times.length if none
     */
    private static int firstIndexAbove(final double[] times, final double value) {
        int low = 0;
        int high = times.length;
        while (low < high) {
            final int mid = (low + high) >>> 1;
            if (times[mid] > value) {
                high = mid;
            } else {
                low = mid + 1;
            }
        }
        return low;
    }

    /**
     * @return the smallest index i with times[i] >= value, or times.length if none
     */
    private static int firstIndexAtOrAbove(final double[] times, final double value) {
        int low = 0;
        int high = times.length;
        while (low < high) {
            final int mid = (low + high) >>> 1;
            if (times[mid] >= value) {
                high = mid;
            } else {
                low = mid + 1;
            }
        }
        return low;
    }

    /**
     * Counting sort into decreasing level, stable within a level.
     */
    private void sortNodeOperationsByReverseLevel() {
        int maxLevel = 0;
        for (DataLikelihoodDelegate.NodeOperation op : nodeOperations) {
            maxLevel = Math.max(maxLevel, op.getLevel());
        }

        // start[l] is where the operations at level l go; the deepest level comes first
        final int[] start = new int[maxLevel + 2];
        for (DataLikelihoodDelegate.NodeOperation op : nodeOperations) {
            ++start[op.getLevel()];
        }
        int position = 0;
        for (int level = maxLevel; level >= 0; --level) {
            final int size = start[level];
            start[level] = position;
            position += size;
        }

        final DataLikelihoodDelegate.NodeOperation[] sorted =
                new DataLikelihoodDelegate.NodeOperation[nodeOperations.size()];
        for (DataLikelihoodDelegate.NodeOperation op : nodeOperations) {
            sorted[start[op.getLevel()]++] = op;
        }

        nodeOperations.clear();
        nodeOperations.addAll(Arrays.asList(sorted));
    }
}
