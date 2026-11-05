/*
 * RootChangingSprMove.java
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

import java.util.ArrayList;
import java.util.BitSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

/**
 * The SPR moves that give the root another node (plan 7 §5.4), Wilson-Balding's two root cases
 * ({@link WilsonBalding#pruneAndRegraft}):
 * <ul>
 * <li><b>above the root</b>: the subtree below i, with its non-root parent P, is pruned (P's other child C joins P's
 * parent Q) and regrafted above the root ρ, so P becomes the root, at h(ρ) + δ with δ the mean length of the root's
 * branches (at least 1% of h(ρ));</li>
 * <li><b>prune at the root</b>: i is a child of the root P; P's other child C becomes the root, and P is regrafted on
 * the branch above j, a node below C whose parent k is older than i, at the middle of (max(h(i), h(j)), h(k)).</li>
 * </ul>
 * TreeModel keeps node numbers with node objects, so heights and moves still refer to nodes by number; the per-branch
 * rate parameters, whose indices skip the root's number, are shifted by their own models when the edit ends
 * (TreeParameterModel.handleRootMove), so every node that keeps a branch keeps its rate. The new branch (ρ's above
 * the root, P's below it) takes the rate that was at the new root's index, which is the rate of the new root's old
 * branch. With the substitution-preserving warm start (§5.6) the changed branches are rescaled as in a clamped SPR.
 *
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public final class RootChangingSprMove implements TimeTreeMove {

    public enum Kind {
        ABOVE_ROOT, PRUNE_AT_ROOT
    }

    private static final double MIN_DELTA = 0.01;

    private final Kind kind;
    private final int i;
    private final int p;
    private final int c;
    private final int q; // P's parent (above the root); -1 (prune at the root)
    private final int j; // the root (above the root); the new place's lower node (prune at the root)
    private final int k; // -1 (above the root); j's parent (prune at the root)

    private RootChangingSprMove(Kind kind, int i, int p, int c, int q, int j, int k) {
        this.kind = kind;
        this.i = i;
        this.p = p;
        this.c = c;
        this.q = q;
        this.j = j;
        this.k = k;
    }

    /**
     * @return the move of the subtree below i (whose parent is not the root) above the root
     */
    public static RootChangingSprMove aboveRoot(TreeModel tree, int i) {
        final NodeRef nodeI = tree.getNode(i);
        final NodeRef nodeP = tree.getParent(nodeI);
        return new RootChangingSprMove(Kind.ABOVE_ROOT, i, nodeP.getNumber(), other(tree, nodeP, nodeI).getNumber(),
                tree.getParent(nodeP).getNumber(), tree.getRoot().getNumber(), -1);
    }

    /**
     * @return the move of the root child i, with the root, to the branch above j (below i's sibling)
     */
    public static RootChangingSprMove pruneAtRoot(TreeModel tree, int i, int j) {
        final NodeRef nodeI = tree.getNode(i);
        final NodeRef root = tree.getRoot();
        return new RootChangingSprMove(Kind.PRUNE_AT_ROOT, i, root.getNumber(), other(tree, root, nodeI).getNumber(),
                -1, j, tree.getParent(tree.getNode(j)).getNumber());
    }

    /**
     * Every above-root move, then every prune-at-root move whose target is at most maxDepth branches below the new
     * root, that the constraints allow, in node-number order.
     */
    public static List<RootChangingSprMove> enumerate(TreeModel tree, int maxDepth, TopologyConstraints constraints) {
        final List<RootChangingSprMove> moves = new ArrayList<RootChangingSprMove>();
        final BitSet[] clades = constraints != null && constraints.size() > 0 ? TopologyConstraints.clades(tree) : null;
        for (int n = 0; n < tree.getNodeCount(); ++n) {
            final NodeRef node = tree.getNode(n);
            if (!tree.isRoot(node) && !tree.isRoot(tree.getParent(node))) {
                add(tree, aboveRoot(tree, n), constraints, clades, moves);
            }
        }
        final NodeRef root = tree.getRoot();
        for (int r = 0; r < 2; ++r) {
            final NodeRef nodeI = tree.getChild(root, r);
            final double hI = tree.getNodeHeight(nodeI);
            final List<NodeRef> level = new ArrayList<NodeRef>();
            level.add(tree.getChild(root, 1 - r));
            // breadth-first below the new root; a branch whose parent is not older than i stops the descent
            for (int depth = 1; depth <= maxDepth && !level.isEmpty(); ++depth) {
                final List<NodeRef> next = new ArrayList<NodeRef>();
                for (NodeRef parent : level) {
                    if (!(tree.getNodeHeight(parent) > hI)) {
                        continue;
                    }
                    for (int ch = 0; ch < tree.getChildCount(parent); ++ch) {
                        final NodeRef nodeJ = tree.getChild(parent, ch);
                        add(tree, pruneAtRoot(tree, nodeI.getNumber(), nodeJ.getNumber()), constraints, clades,
                                moves);
                        next.add(nodeJ);
                    }
                }
                level.clear();
                level.addAll(next);
            }
        }
        return moves;
    }

    private static void add(TreeModel tree, RootChangingSprMove move, TopologyConstraints constraints,
                            BitSet[] clades, List<RootChangingSprMove> moves) {
        if (clades == null || !constraints.violates(tree, move, clades)) {
            moves.add(move);
        }
    }

    public Kind getKind() {
        return kind;
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

    public int getK() {
        return k;
    }

    @Override
    public boolean isValid(TreeModel tree) {
        final NodeRef nodeI = tree.getNode(i);
        final NodeRef nodeP = tree.getNode(p);
        if (tree.isRoot(nodeI) || tree.getParent(nodeI) != nodeP || other(tree, nodeP, nodeI).getNumber() != c) {
            return false;
        }
        if (kind == Kind.ABOVE_ROOT) {
            return !tree.isRoot(nodeP) && tree.getParent(nodeP).getNumber() == q && tree.getRoot().getNumber() == j;
        }
        final NodeRef nodeJ = tree.getNode(j);
        if (!tree.isRoot(nodeP) || tree.isRoot(nodeJ) || tree.getParent(nodeJ).getNumber() != k || j == c) {
            return false;
        }
        for (NodeRef x = nodeJ; x != null; x = tree.getParent(x)) {
            if (x.getNumber() == c) {
                return tree.getNodeHeight(tree.getNode(k)) > tree.getNodeHeight(nodeI);
            }
        }
        return false;
    }

    // a root move always sets P's height
    @Override
    public boolean needsRepair(TreeModel tree) {
        return true;
    }

    /**
     * @return P's height after the move
     */
    public double newHeight(TreeModel tree) {
        if (kind == Kind.ABOVE_ROOT) {
            final NodeRef root = tree.getRoot();
            final double hRoot = tree.getNodeHeight(root);
            final double mean = 0.5 * (2.0 * hRoot - tree.getNodeHeight(tree.getChild(root, 0))
                    - tree.getNodeHeight(tree.getChild(root, 1)));
            return hRoot + Math.max(mean, MIN_DELTA * hRoot);
        }
        final double lo = Math.max(tree.getNodeHeight(tree.getNode(i)), tree.getNodeHeight(tree.getNode(j)));
        return lo + 0.5 * (tree.getNodeHeight(tree.getNode(k)) - lo);
    }

    @Override
    public void apply(TreeModel tree, BranchRateWarmStart warmStart) {
        if (!isValid(tree)) {
            throw new IllegalStateException("The tree no longer has the nodes of " + describe());
        }
        final NodeRef nodeI = tree.getNode(i);
        final NodeRef nodeP = tree.getNode(p);
        final NodeRef nodeC = tree.getNode(c);
        final double hP = tree.getNodeHeight(nodeP);
        final double hI = tree.getNodeHeight(nodeI);
        final double newP = newHeight(tree);

        // the new rates (plan 7 §5.6), from the old ones, before the edit
        double rateC = 0, rateP = 0, rateRoot = 0;
        final boolean warm = warmStart != null && warmStart.isActive();
        final double rateI = warm ? warmStart.preserve(warmStart.getRate(nodeI), hP - hI, newP - hI) : 0;
        if (warm) {
            if (kind == Kind.ABOVE_ROOT) {
                final double hC = tree.getNodeHeight(nodeC);
                final double hQ = tree.getNodeHeight(tree.getNode(q));
                final double rC = warmStart.getRate(nodeC);
                // C's branch now runs to Q, as in a clamped SPR; the old root's new branch takes P's old rate
                final double merged = (rC * (hP - hC) + warmStart.getRate(nodeP) * (hQ - hP)) / (hQ - hC);
                rateC = warmStart.preserve(rC, merged, rC);
                rateRoot = warmStart.getRate(nodeP);
            } else {
                // P's branch is the upper part of j's old branch
                rateP = warmStart.getRate(tree.getNode(j));
            }
        }

        tree.beginTreeEdit();
        WilsonBalding.pruneAndRegraft(tree, nodeI, tree.getNode(j));
        tree.setNodeHeight(nodeP, newP);
        tree.endTreeEdit();

        if (warm) {
            warmStart.setRate(tree.getNode(i), rateI);
            if (kind == Kind.ABOVE_ROOT) {
                warmStart.setRate(tree.getNode(c), rateC);
                warmStart.setRate(tree.getNode(j), rateRoot);
            } else {
                warmStart.setRate(tree.getNode(p), rateP);
            }
        }
    }

    @Override
    public void applyTopology(TreeModel tree) {
        WilsonBalding.pruneAndRegraft(tree, tree.getNode(i), tree.getNode(j));
    }

    /**
     * @return above the root: P, the old root and Q; at the root: P, k and C (the new root)
     */
    @Override
    public int[] localHeights(TreeModel tree) {
        return kind == Kind.ABOVE_ROOT ? toArray(p, j, q) : toArray(p, k, c);
    }

    /**
     * @return the branches whose length changed: above the root, those above i, C and the old root; at the root, those
     * above i, j and P
     */
    @Override
    public int[] localRates(TreeModel tree) {
        return kind == Kind.ABOVE_ROOT ? toArray(i, c, j) : toArray(i, j, p);
    }

    @Override
    public int[] touchedNodes() {
        return kind == Kind.ABOVE_ROOT ? toArray(i, p, c, q, j) : toArray(i, p, c, j, k);
    }

    private static int[] toArray(int... nodes) {
        Set<Integer> set = new LinkedHashSet<Integer>();
        for (int node : nodes) {
            set.add(node);
        }
        int[] result = new int[set.size()];
        int n = 0;
        for (int x : set) {
            result[n++] = x;
        }
        return result;
    }

    private static NodeRef other(TreeModel tree, NodeRef parent, NodeRef child) {
        return tree.getChild(parent, 0) == child ? tree.getChild(parent, 1) : tree.getChild(parent, 0);
    }

    @Override
    public String describe() {
        return kind == Kind.ABOVE_ROOT ? "rootSPR(i=" + i + ", P=" + p + " -> above the root " + j + ")" :
                "rootSPR(i=" + i + ", root P=" + p + " -> above j=" + j + ", C=" + c + " the root)";
    }

    @Override
    public boolean equals(Object o) {
        if (!(o instanceof RootChangingSprMove)) {
            return false;
        }
        RootChangingSprMove m = (RootChangingSprMove) o;
        return kind == m.kind && i == m.i && p == m.p && c == m.c && q == m.q && j == m.j && k == m.k;
    }

    @Override
    public int hashCode() {
        return (((((kind.ordinal() * 31 + i) * 31 + p) * 31 + c) * 31 + q) * 31 + j) * 31 + k;
    }

    @Override
    public String toString() {
        return describe();
    }
}
