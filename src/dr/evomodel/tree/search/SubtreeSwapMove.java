/*
 * SubtreeSwapMove.java
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

package dr.evomodel.tree.search;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.operators.AbstractTreeOperator;
import dr.evomodel.tree.TreeModel;

import java.util.ArrayList;
import java.util.List;

/**
 * A subtree swap at fixed heights, BEAST's wide exchange (ExchangeOperator.wide, GibbsSubtreeSwap): the subtrees below
 * i and j change parents, with no height changed. It reaches trees two SPRs away in one move. The edit is
 * AbstractTreeOperator.exchange; the feasibility test and the node distance are those of ExchangeOperator.wide and
 * GibbsSubtreeSwap (which keep private copies of the distance), here until they are shared. Nodes are referred to by
 * number.
 *
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public final class SubtreeSwapMove implements TimeTreeMove {

    private final int i;
    private final int iP;
    private final int j;
    private final int jP;

    public SubtreeSwapMove(int i, int iP, int j, int jP) {
        this.i = i;
        this.iP = iP;
        this.j = j;
        this.jP = jP;
    }

    /**
     * @return the swaps whose parents are between minDistance and maxDistance nodes apart (as GibbsSubtreeSwap's
     * window; distance 1 is a narrow exchange, an NNI), each pair once, by node number
     */
    public static List<SubtreeSwapMove> enumerate(TreeModel tree, int minDistance, int maxDistance,
                                                  TopologyConstraints constraints) {
        final java.util.BitSet[] clades = constraints != null && constraints.size() > 0 ?
                TopologyConstraints.clades(tree) : null;
        List<SubtreeSwapMove> moves = new ArrayList<SubtreeSwapMove>();
        for (int a = 0; a < tree.getNodeCount(); ++a) {
            final NodeRef nodeA = tree.getNode(a);
            if (tree.isRoot(nodeA)) {
                continue;
            }
            for (int b = a + 1; b < tree.getNodeCount(); ++b) {
                final NodeRef nodeB = tree.getNode(b);
                if (tree.isRoot(nodeB) || !canExchange(tree, nodeA, nodeB)) {
                    continue;
                }
                final int d = nodeDistance(tree, tree.getParent(nodeA), tree.getParent(nodeB));
                if (d < minDistance || d > maxDistance) {
                    continue;
                }
                SubtreeSwapMove move = new SubtreeSwapMove(a, tree.getParent(nodeA).getNumber(), b,
                        tree.getParent(nodeB).getNumber());
                if (clades == null || !constraints.violates(tree, move, clades)) {
                    moves.add(move);
                }
            }
        }
        return moves;
    }

    /**
     * @return whether the subtrees below non-root nodes i and j can be exchanged at fixed heights, as in
     * ExchangeOperator.wide: different parents, neither the other's parent, each younger than the other's parent
     * (which also rules out one inside the other)
     */
    static boolean canExchange(Tree tree, NodeRef i, NodeRef j) {
        final NodeRef iP = tree.getParent(i);
        final NodeRef jP = tree.getParent(j);
        return iP != jP && i != jP && j != iP
                && tree.getNodeHeight(j) < tree.getNodeHeight(iP)
                && tree.getNodeHeight(i) < tree.getNodeHeight(jP);
    }

    /**
     * @return the nodes between i and j, as getNodeDistance in ExchangeOperator and GibbsSubtreeSwap: each step moves
     * the younger of the two up to its parent
     */
    static int nodeDistance(Tree tree, NodeRef i, NodeRef j) {
        int count = 0;
        double heightI = tree.getNodeHeight(i);
        double heightJ = tree.getNodeHeight(j);
        while (i != j) {
            count++;
            if (heightI < heightJ) {
                i = tree.getParent(i);
                heightI = tree.getNodeHeight(i);
            } else {
                j = tree.getParent(j);
                heightJ = tree.getNodeHeight(j);
            }
        }
        return count;
    }

    public int getI() {
        return i;
    }

    public int getIP() {
        return iP;
    }

    public int getJ() {
        return j;
    }

    public int getJP() {
        return jP;
    }

    @Override
    public boolean isValid(TreeModel tree) {
        final NodeRef nodeI = tree.getNode(i);
        final NodeRef nodeJ = tree.getNode(j);
        return !tree.isRoot(nodeI) && !tree.isRoot(nodeJ) && tree.getParent(nodeI).getNumber() == iP
                && tree.getParent(nodeJ).getNumber() == jP && canExchange(tree, nodeI, nodeJ);
    }

    @Override
    public boolean needsRepair(TreeModel tree) {
        return false;
    }

    @Override
    public void apply(TreeModel tree, BranchRateWarmStart warmStart) {
        if (!isValid(tree)) {
            throw new IllegalStateException("The tree no longer has the nodes of " + describe());
        }
        final NodeRef nodeI = tree.getNode(i);
        final NodeRef nodeJ = tree.getNode(j);
        final double hI = tree.getNodeHeight(nodeI);
        final double hJ = tree.getNodeHeight(nodeJ);
        final double hIP = tree.getNodeHeight(tree.getNode(iP));
        final double hJP = tree.getNodeHeight(tree.getNode(jP));
        double rateI = 0, rateJ = 0;
        final boolean warm = warmStart != null && warmStart.isActive();
        if (warm) {
            rateI = warmStart.preserve(warmStart.getRate(nodeI), hIP - hI, hJP - hI);
            rateJ = warmStart.preserve(warmStart.getRate(nodeJ), hJP - hJ, hIP - hJ);
        }
        tree.beginTreeEdit();
        applyTopology(tree);
        tree.endTreeEdit();
        if (warm) {
            warmStart.setRate(tree.getNode(i), rateI);
            warmStart.setRate(tree.getNode(j), rateJ);
        }
    }

    @Override
    public void applyTopology(TreeModel tree) {
        final NodeRef nodeI = tree.getNode(i);
        final NodeRef nodeJ = tree.getNode(j);
        AbstractTreeOperator.exchange(tree, nodeI, nodeJ, tree.getParent(nodeI), tree.getParent(nodeJ));
    }

    /**
     * @return the two parents, whose children changed
     */
    @Override
    public int[] localHeights(TreeModel tree) {
        return new int[]{iP, jP};
    }

    /**
     * @return the branches above i and j, which changed length, and above their parents (but the root's)
     */
    @Override
    public int[] localRates(TreeModel tree) {
        List<Integer> nodes = new ArrayList<Integer>();
        nodes.add(i);
        nodes.add(j);
        for (int parent : new int[]{iP, jP}) {
            if (!tree.isRoot(tree.getNode(parent))) {
                nodes.add(parent);
            }
        }
        int[] result = new int[nodes.size()];
        for (int k = 0; k < result.length; ++k) {
            result[k] = nodes.get(k);
        }
        return result;
    }

    @Override
    public int[] touchedNodes() {
        return new int[]{i, iP, j, jP};
    }

    @Override
    public String describe() {
        return "swap(" + i + " <-> " + j + ")";
    }

    @Override
    public boolean equals(Object o) {
        if (!(o instanceof SubtreeSwapMove)) {
            return false;
        }
        SubtreeSwapMove m = (SubtreeSwapMove) o;
        return i == m.i && iP == m.iP && j == m.j && jP == m.jP;
    }

    @Override
    public int hashCode() {
        return ((i * 31 + iP) * 31 + j) * 31 + jP;
    }

    @Override
    public String toString() {
        return describe();
    }
}
