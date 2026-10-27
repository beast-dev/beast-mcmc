/*
 * RootedNniMove.java
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

/**
 * A rooted nearest-neighbour interchange with height repair (plan 7 §5.1). At a non-root internal node P with parent
 * G, children a and b, and sibling U, swap a and U: a moves under G and U under P. The result is feasible as long as
 * U is younger than P; otherwise P is moved to the middle of its new interval, m + (h(G) - m) / 2 with
 * m = max(h(b), h(U)), which is never empty. So every rooted NNI can be made feasible, and the moves at all non-root
 * internal nodes number 2 (n - 2).
 * <p>
 * The rate warm start keeps each rewritten branch's expected substitutions, rate times time (§5.6): a's and U's
 * branches change length, and with a repair also b's and P's.
 *
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public final class RootedNniMove implements TimeTreeMove {

    private final int p;
    private final int g;
    private final int a;
    private final int u;

    public RootedNniMove(int p, int g, int a, int u) {
        this.p = p;
        this.g = g;
        this.a = a;
        this.u = u;
    }

    public int getP() {
        return p;
    }

    public int getG() {
        return g;
    }

    public int getA() {
        return a;
    }

    public int getU() {
        return u;
    }

    /**
     * @return b, P's child other than a, before the move
     */
    public int getB(TreeModel tree) {
        return other(tree).getNumber();
    }

    /**
     * @return the quartet's five nodes before the move: P, G, a, b and U
     */
    public int[] quartet(TreeModel tree) {
        return new int[]{p, g, a, getB(tree), u};
    }

    /**
     * @return what a scored move's value depends on besides its own nodes, before the move: b and G's parent (-1 at
     * the root). A move kept from an earlier round is re-scored only while these are unchanged
     */
    public int[] neighbourhood(TreeModel tree) {
        final NodeRef nodeG = tree.getNode(g);
        return new int[]{getB(tree), tree.isRoot(nodeG) ? -1 : tree.getParent(nodeG).getNumber()};
    }

    @Override
    public boolean equals(Object o) {
        if (!(o instanceof RootedNniMove)) {
            return false;
        }
        RootedNniMove m = (RootedNniMove) o;
        return p == m.p && g == m.g && a == m.a && u == m.u;
    }

    @Override
    public int hashCode() {
        return ((p * 31 + g) * 31 + a) * 31 + u;
    }

    @Override
    public boolean isValid(TreeModel tree) {
        final NodeRef nodeP = tree.getNode(p);
        final NodeRef nodeG = tree.getNode(g);
        return !tree.isRoot(nodeP) && tree.getParent(nodeP) == nodeG && tree.getParent(tree.getNode(a)) == nodeP
                && tree.getParent(tree.getNode(u)) == nodeG && u != p;
    }

    @Override
    public boolean needsRepair(TreeModel tree) {
        return !(tree.getNodeHeight(tree.getNode(u)) < tree.getNodeHeight(tree.getNode(p)));
    }

    /**
     * @return P's height after the move: unchanged when U is younger than P, else the middle of its new interval
     */
    public double repairedHeight(TreeModel tree) {
        final double hP = tree.getNodeHeight(tree.getNode(p));
        if (!needsRepair(tree)) {
            return hP;
        }
        final double m = Math.max(tree.getNodeHeight(other(tree)), tree.getNodeHeight(tree.getNode(u)));
        return m + 0.5 * (tree.getNodeHeight(tree.getNode(g)) - m);
    }

    @Override
    public void apply(TreeModel tree, BranchRateWarmStart warmStart) {
        if (!isValid(tree)) {
            throw new IllegalStateException("The tree no longer has the nodes of " + describe());
        }
        final NodeRef nodeP = tree.getNode(p);
        final NodeRef nodeA = tree.getNode(a);
        final NodeRef nodeU = tree.getNode(u);
        final NodeRef nodeB = other(tree);
        final int b = nodeB.getNumber();
        final double hP = tree.getNodeHeight(nodeP);
        final double hG = tree.getNodeHeight(tree.getNode(g));
        final double hA = tree.getNodeHeight(nodeA);
        final double hU = tree.getNodeHeight(nodeU);
        final double hB = tree.getNodeHeight(nodeB);
        final boolean repair = needsRepair(tree);
        final double newP = repairedHeight(tree);

        // the new rates, from the old ones, before the edit
        double rateA = 0, rateU = 0, rateB = 0, rateP = 0;
        final boolean warm = warmStart != null && warmStart.isActive();
        if (warm) {
            rateA = warmStart.preserve(warmStart.getRate(nodeA), hP - hA, hG - hA);
            rateU = warmStart.preserve(warmStart.getRate(nodeU), hG - hU, newP - hU);
            if (repair) {
                rateB = warmStart.preserve(warmStart.getRate(nodeB), hP - hB, newP - hB);
                rateP = warmStart.preserve(warmStart.getRate(nodeP), hG - hP, hG - newP);
            }
        }

        tree.beginTreeEdit();
        applyTopology(tree);
        if (repair) {
            tree.setNodeHeight(nodeP, newP);
        }
        tree.endTreeEdit();

        if (warm) {
            warmStart.setRate(tree.getNode(a), rateA);
            warmStart.setRate(tree.getNode(u), rateU);
            if (repair) {
                warmStart.setRate(tree.getNode(b), rateB);
                warmStart.setRate(tree.getNode(p), rateP);
            }
        }
    }

    // b after the move, when a has left P: P's child other than U
    private int otherAfterMove(TreeModel tree) {
        final NodeRef nodeP = tree.getNode(p);
        final NodeRef first = tree.getChild(nodeP, 0);
        return first.getNumber() == u ? tree.getChild(nodeP, 1).getNumber() : first.getNumber();
    }

    @Override
    public void applyTopology(TreeModel tree) {
        TreeEdits.exchange(tree, tree.getNode(a), tree.getNode(u));
    }

    // b: P's child other than a, before the move
    private NodeRef other(TreeModel tree) {
        final NodeRef nodeP = tree.getNode(p);
        final NodeRef first = tree.getChild(nodeP, 0);
        return first.getNumber() == a ? tree.getChild(nodeP, 1) : first;
    }

    /**
     * @return P and G: the heights of the quartet's two internal nodes
     */
    @Override
    public int[] localHeights(TreeModel tree) {
        return new int[]{p, g};
    }

    /**
     * @return the five branches around the quartet after the move: a, U, b, P and G's own (unless G is the root)
     */
    @Override
    public int[] localRates(TreeModel tree) {
        final int b = tree.getParent(tree.getNode(a)).getNumber() == g ? otherAfterMove(tree) : other(tree).getNumber();
        return tree.isRoot(tree.getNode(g)) ? new int[]{a, u, b, p} : new int[]{a, u, b, p, g};
    }

    @Override
    public int[] touchedNodes() {
        return new int[]{p, g, a, u};
    }

    @Override
    public String describe() {
        return "NNI(P=" + p + ", G=" + g + ", a=" + a + " <-> U=" + u + ")";
    }

    @Override
    public String toString() {
        return describe();
    }
}
