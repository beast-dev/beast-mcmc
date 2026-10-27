/*
 * SprNeighbourhood.java
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
import dr.evomodel.tree.TreeModel;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

/**
 * The clamped SPR moves of a pruned subtree (plan 7 §5.3): the branches of the pruned tree within a window of
 * distances from the pruned subtree's old place, found breadth-first, scored as they are found. As in RAxML's lazy
 * SPR, a target's neighbours are visited only while its score stays within a cutoff of the incumbent's; branches whose
 * parent is not older than the pruned subtree's root cannot take it, and neither can any branch below them, so from
 * those the search only goes up. Distances count branches in the pruned tree, where C's branch now runs to Q; the two
 * branches below the root are one apart.
 *
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public final class SprNeighbourhood {

    /**
     * Scores a target as it is found.
     */
    public interface Visitor {
        /**
         * @return the move's score, compared with the threshold to decide whether to look beyond its target
         */
        double visit(ClampedSprMove move, int distance);
    }

    private SprNeighbourhood() {
    }

    /**
     * @return the nodes whose subtree can be pruned (its parent is not the root), in priority order: increasing
     * relative branch time (h(p) - h(i)) / h(p), ties by node number
     */
    public static List<Integer> pruneOrder(final TreeModel tree) {
        List<Integer> nodes = new ArrayList<Integer>();
        final double[] priority = new double[tree.getNodeCount()];
        for (int n = 0; n < tree.getNodeCount(); ++n) {
            NodeRef node = tree.getNode(n);
            if (!tree.isRoot(node) && !tree.isRoot(tree.getParent(node))) {
                nodes.add(n);
                final double parent = tree.getNodeHeight(tree.getParent(node));
                priority[n] = parent > 0.0 ? (parent - tree.getNodeHeight(node)) / parent : 0.0;
            }
        }
        Collections.sort(nodes, new Comparator<Integer>() {
            public int compare(Integer x, Integer y) {
                int c = Double.compare(priority[x], priority[y]);
                return c != 0 ? c : Integer.compare(x, y);
            }
        });
        return nodes;
    }

    /**
     * @return every target of the subtree below i within [rMin, rMax] that the constraints allow, without scores
     */
    public static List<ClampedSprMove> enumerate(TreeModel tree, int i, int rMin, int rMax,
                                                 TopologyConstraints constraints) {
        final List<ClampedSprMove> moves = new ArrayList<ClampedSprMove>();
        search(tree, i, rMin, rMax, constraints, new Visitor() {
            public double visit(ClampedSprMove move, int distance) {
                moves.add(move);
                return Double.POSITIVE_INFINITY;
            }
        }, Double.NEGATIVE_INFINITY);
        return moves;
    }

    /**
     * Visits the targets of the subtree below i (whose parent must not be the root) at distances in [rMin, rMax],
     * breadth-first from its old place. A visited target's neighbours are looked at only if its score is at least the
     * threshold; branches closer than rMin, or that the constraints forbid, are passed through without a score. The
     * visitor may score moves on the tree (store, move, restore): restores swap TreeModel's node objects, so the
     * search holds node numbers only.
     */
    public static void search(TreeModel tree, int i, int rMin, int rMax, TopologyConstraints constraints,
                              Visitor visitor, double threshold) {
        final Pruned pruned = new Pruned(tree, i);
        final double hI = tree.getNodeHeight(tree.getNode(i));
        final BitSet[] clades = constraints != null && constraints.size() > 0 ? TopologyConstraints.clades(tree) : null;

        final BitSet visited = new BitSet(tree.getNodeCount());
        visited.set(pruned.c);
        final List<int[]> queue = new ArrayList<int[]>(); // node, distance, 1 if only up from it
        queue.add(new int[]{pruned.c, 0, 0});
        for (int head = 0; head < queue.size(); ++head) {
            final int x = queue.get(head)[0];
            final int d = queue.get(head)[1];
            if (d >= rMax) {
                continue;
            }
            // the branches adjacent to x's: below x (unless only going up), x's siblings, and the one above them
            final List<Integer> next = new ArrayList<Integer>();
            if (queue.get(head)[2] == 0) {
                next.addAll(pruned.children(x));
            }
            final int u = pruned.parent(x);
            for (int sibling : pruned.children(u)) {
                if (sibling != x) {
                    next.add(sibling);
                }
            }
            if (pruned.parent(u) >= 0) {
                next.add(u);
            }
            for (int y : next) {
                if (visited.get(y)) {
                    continue;
                }
                visited.set(y);
                // a branch no older above than i cannot take it, nor can any below it: only up from there
                final boolean feasible = tree.getNodeHeight(tree.getNode(pruned.parent(y))) > hI;
                if (feasible && d + 1 >= rMin) {
                    ClampedSprMove move = new ClampedSprMove(i, pruned.p, pruned.c, pruned.q, y,
                            tree.getParent(tree.getNode(y)).getNumber());
                    if ((clades == null || !constraints.violates(tree, move, clades))
                            && !(visitor.visit(move, d + 1) >= threshold)) {
                        continue; // scored below the cutoff: nothing beyond it
                    }
                }
                queue.add(new int[]{y, d + 1, feasible ? 0 : 1});
            }
        }
    }

    // the tree with P removed, by node number: C's parent is Q, and Q's children are C and P's sibling
    private static final class Pruned {
        private final TreeModel tree;
        final int p;
        final int c;
        final int q;

        Pruned(TreeModel tree, int i) {
            this.tree = tree;
            final NodeRef nodeI = tree.getNode(i);
            final NodeRef nodeP = tree.getParent(nodeI);
            this.p = nodeP.getNumber();
            this.c = (tree.getChild(nodeP, 0) == nodeI ? tree.getChild(nodeP, 1) : tree.getChild(nodeP, 0)).getNumber();
            this.q = tree.getParent(nodeP).getNumber();
        }

        int parent(int x) {
            if (x == c) {
                return q;
            }
            final NodeRef node = tree.getNode(x);
            return tree.isRoot(node) ? -1 : tree.getParent(node).getNumber();
        }

        List<Integer> children(int x) {
            final NodeRef node = tree.getNode(x);
            List<Integer> children = new ArrayList<Integer>();
            for (int k = 0; k < tree.getChildCount(node); ++k) {
                final int child = tree.getChild(node, k).getNumber();
                children.add(child == p ? c : child);
            }
            return children;
        }
    }
}
