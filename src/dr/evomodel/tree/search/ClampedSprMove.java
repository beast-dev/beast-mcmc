/*
 * ClampedSprMove.java
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
import dr.evomodel.operators.WilsonBalding;
import dr.evomodel.tree.TreeModel;

import java.util.LinkedHashSet;
import java.util.Set;

/**
 * Clamped-height SPR (plan 7 §5.3), a generalised fixed-height SPR: the subtree below i, with its parent P (not the
 * root), is pruned (P's other child C joins P's parent Q) and regrafted on the branch above j, whose parent is pj.
 * P keeps its height when it lies strictly between max(h(i), h(j)) and h(pj), as in FixedHeightSubtreePruneRegraft;
 * otherwise it moves to the middle of that interval, as in Wilson-Balding with the draw replaced by the midpoint. The
 * edit itself is WilsonBalding.pruneAndRegraft. Nodes are referred to by number.
 *
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public final class ClampedSprMove implements TimeTreeMove {

    private final int i;
    private final int p;
    private final int c;
    private final int q;
    private final int j;
    private final int pj;

    /**
     * @param i  the root of the pruned subtree
     * @param p  its parent, not the root
     * @param c  p's other child
     * @param q  p's parent
     * @param j  the node above which the subtree is regrafted, outside subtree(i), not c
     * @param pj j's parent, not P
     */
    public ClampedSprMove(int i, int p, int c, int q, int j, int pj) {
        this.i = i;
        this.p = p;
        this.c = c;
        this.q = q;
        this.j = j;
        this.pj = pj;
    }

    /**
     * @return the move of the subtree below i to the branch above j, from the tree as it is
     */
    public static ClampedSprMove of(TreeModel tree, int i, int j) {
        final NodeRef nodeI = tree.getNode(i);
        final NodeRef nodeP = tree.getParent(nodeI);
        final NodeRef nodeC = tree.getChild(nodeP, 0) == nodeI ? tree.getChild(nodeP, 1) : tree.getChild(nodeP, 0);
        return new ClampedSprMove(i, nodeP.getNumber(), nodeC.getNumber(), tree.getParent(nodeP).getNumber(), j,
                tree.getParent(tree.getNode(j)).getNumber());
    }

    public int getI() {
        return i;
    }

    public int getP() {
        return p;
    }

    public int getC() {
        return c;
    }

    public int getQ() {
        return q;
    }

    public int getJ() {
        return j;
    }

    public int getPj() {
        return pj;
    }

    @Override
    public boolean isValid(TreeModel tree) {
        final NodeRef nodeI = tree.getNode(i);
        final NodeRef nodeP = tree.getNode(p);
        if (tree.isRoot(nodeI) || tree.getParent(nodeI) != nodeP || tree.isRoot(nodeP)
                || tree.getParent(nodeP).getNumber() != q || other(tree, nodeP, nodeI).getNumber() != c) {
            return false;
        }
        final NodeRef nodeJ = tree.getNode(j);
        return !tree.isRoot(nodeJ) && tree.getParent(nodeJ).getNumber() == pj && j != c && j != p && pj != p
                && !isBelow(tree, nodeJ, nodeI);
    }

    private static boolean isBelow(TreeModel tree, NodeRef node, NodeRef ancestor) {
        for (NodeRef x = node; x != null; x = tree.getParent(x)) {
            if (x == ancestor) {
                return true;
            }
        }
        return false;
    }

    private static NodeRef other(TreeModel tree, NodeRef parent, NodeRef child) {
        return tree.getChild(parent, 0) == child ? tree.getChild(parent, 1) : tree.getChild(parent, 0);
    }

    /**
     * @return whether P's height must change: it is not strictly between max(h(i), h(j)) and h(pj)
     */
    @Override
    public boolean needsRepair(TreeModel tree) {
        final double h = tree.getNodeHeight(tree.getNode(p));
        return !(lower(tree) < h && h < tree.getNodeHeight(tree.getNode(pj)));
    }

    private double lower(TreeModel tree) {
        return Math.max(tree.getNodeHeight(tree.getNode(i)), tree.getNodeHeight(tree.getNode(j)));
    }

    /**
     * @return P's height after the move: kept when feasible, else the middle of its new interval
     */
    public double newHeight(TreeModel tree) {
        if (!needsRepair(tree)) {
            return tree.getNodeHeight(tree.getNode(p));
        }
        final double lo = lower(tree);
        return lo + 0.5 * (tree.getNodeHeight(tree.getNode(pj)) - lo);
    }

    /**
     * @return whether the new branch above P exists at all: pj older than i
     */
    public boolean isFeasible(TreeModel tree) {
        return tree.getNodeHeight(tree.getNode(pj)) > tree.getNodeHeight(tree.getNode(i));
    }

    @Override
    public void apply(TreeModel tree, BranchRateWarmStart warmStart) {
        if (!isValid(tree)) {
            throw new IllegalStateException("The tree no longer has the nodes of " + describe());
        }
        final NodeRef nodeI = tree.getNode(i);
        final NodeRef nodeP = tree.getNode(p);
        final NodeRef nodeC = tree.getNode(c);
        final NodeRef nodeJ = tree.getNode(j);
        final double hP = tree.getNodeHeight(nodeP);
        final double hC = tree.getNodeHeight(nodeC);
        final double hQ = tree.getNodeHeight(tree.getNode(q));
        final double hI = tree.getNodeHeight(nodeI);
        final boolean repair = needsRepair(tree);
        final double newP = newHeight(tree);

        // the new rates (plan 7 §5.6), from the old ones, before the edit
        double rateC = 0, rateP = 0, rateI = 0;
        final boolean warm = warmStart != null && warmStart.isActive();
        if (warm) {
            final double rC = warmStart.getRate(nodeC);
            // C's branch now runs to Q: the substitutions of C's and P's old branches over its new length
            final double merged = (rC * (hP - hC) + warmStart.getRate(nodeP) * (hQ - hP)) / (hQ - hC);
            rateC = warmStart.preserve(rC, merged, rC);
            // P's branch is the upper part of j's old branch
            rateP = warmStart.getRate(nodeJ);
            if (repair) {
                rateI = warmStart.preserve(warmStart.getRate(nodeI), hP - hI, newP - hI);
            }
        }

        tree.beginTreeEdit();
        WilsonBalding.pruneAndRegraft(tree, nodeI, nodeJ);
        if (repair) {
            tree.setNodeHeight(nodeP, newP);
        }
        tree.endTreeEdit();

        if (warm) {
            warmStart.setRate(tree.getNode(c), rateC);
            warmStart.setRate(tree.getNode(p), rateP);
            if (repair) {
                warmStart.setRate(tree.getNode(i), rateI);
            }
        }
    }

    @Override
    public void applyTopology(TreeModel tree) {
        WilsonBalding.pruneAndRegraft(tree, tree.getNode(i), tree.getNode(j));
    }

    /**
     * @return P, Q and pj: P's new height and the heights it lies between (the time-tree triplet)
     */
    @Override
    public int[] localHeights(TreeModel tree) {
        Set<Integer> nodes = new LinkedHashSet<Integer>();
        nodes.add(p);
        nodes.add(q);
        nodes.add(pj);
        return toArray(nodes);
    }

    /**
     * @return the branches above C, i, j and P
     */
    @Override
    public int[] localRates(TreeModel tree) {
        return new int[]{c, i, j, p};
    }

    @Override
    public int[] touchedNodes() {
        Set<Integer> nodes = new LinkedHashSet<Integer>();
        for (int node : new int[]{i, p, c, q, j, pj}) {
            nodes.add(node);
        }
        return toArray(nodes);
    }

    private static int[] toArray(Set<Integer> set) {
        int[] result = new int[set.size()];
        int k = 0;
        for (int x : set) {
            result[k++] = x;
        }
        return result;
    }

    @Override
    public String describe() {
        return "SPR(i=" + i + ", P=" + p + " -> above j=" + j + ")";
    }

    @Override
    public boolean equals(Object o) {
        if (!(o instanceof ClampedSprMove)) {
            return false;
        }
        ClampedSprMove m = (ClampedSprMove) o;
        return i == m.i && p == m.p && c == m.c && q == m.q && j == m.j && pj == m.pj;
    }

    @Override
    public int hashCode() {
        return ((((i * 31 + p) * 31 + c) * 31 + q) * 31 + j) * 31 + pj;
    }

    @Override
    public String toString() {
        return describe();
    }
}
