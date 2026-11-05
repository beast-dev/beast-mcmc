/*
 * TopologyConstraints.java
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
import dr.evolution.tree.TreeUtils;
import dr.evomodel.tree.MonophylyStatistic;
import dr.inference.model.BooleanLikelihood;
import dr.inference.model.BooleanStatistic;
import dr.inference.model.CompoundLikelihood;
import dr.inference.model.Likelihood;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collections;
import java.util.List;
import java.util.Set;

/**
 * Taxon sets that must stay clades during the search (plan 7 §5.7), as tip bitsets: the monophyletic sets of the
 * objective (BEAUti's monophyly statistics, scored by a boolean likelihood) and any set given to the search, such as
 * the taxa a tree prior excludes. Taxon sets that only carry a tMRCA prior are not constraints.
 * <p>
 * A rooted NNI that swaps a child a of P with P's sibling U removes exactly one clade, clade(P), and creates
 * b + U; every other clade stays. So it breaks a constraint S exactly when S = clade(P), which moving a whole
 * constrained clade as a unit never does.
 * <p>
 * Monophyly statistics with taxa that may be ignored, or inverted ones, are not filtered here; the objective still
 * scores their violations as impossible.
 *
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public final class TopologyConstraints {

    private final List<Set<String>> sets = new ArrayList<Set<String>>();
    private final List<BitSet> tips = new ArrayList<BitSet>();
    private int unfiltered = 0;

    /**
     * The constraints of the objective's monophyly statistics, plus the given sets.
     */
    public TopologyConstraints(Tree tree, Likelihood objective, List<Set<String>> extra) {
        if (objective != null) {
            collect(tree, objective);
        }
        if (extra != null) {
            for (Set<String> set : extra) {
                add(tree, set);
            }
        }
    }

    private void collect(Tree tree, Likelihood likelihood) {
        if (likelihood instanceof CompoundLikelihood) {
            for (Likelihood term : ((CompoundLikelihood) likelihood).getLikelihoods()) {
                collect(tree, term);
            }
        } else if (likelihood instanceof BooleanLikelihood) {
            for (BooleanStatistic statistic : ((BooleanLikelihood) likelihood).getData()) {
                if (statistic instanceof MonophylyStatistic) {
                    MonophylyStatistic monophyly = (MonophylyStatistic) statistic;
                    if (monophyly.getTree() == tree && !monophyly.isInverse()
                            && monophyly.getIgnoreLeafSet().isEmpty()) {
                        add(tree, monophyly.getLeafSet());
                    } else {
                        unfiltered++;
                    }
                } else {
                    unfiltered++;
                }
            }
        }
    }

    private void add(Tree tree, Set<String> set) {
        if (set.size() < 2 || set.size() >= tree.getExternalNodeCount()) {
            return; // always a clade
        }
        BitSet bits = new BitSet(tree.getExternalNodeCount());
        for (String id : set) {
            int index = tree.getTaxonIndex(id);
            if (index < 0) {
                throw new IllegalArgumentException("Taxon " + id + " of a constraint is not in the tree");
            }
            bits.set(tipNumber(tree, index));
        }
        if (!tips.contains(bits)) {
            sets.add(Collections.unmodifiableSet(set));
            tips.add(bits);
        }
    }

    // the node number of the taxon's tip
    private static int tipNumber(Tree tree, int taxonIndex) {
        for (int i = 0; i < tree.getExternalNodeCount(); ++i) {
            NodeRef tip = tree.getExternalNode(i);
            if (tree.getNodeTaxon(tip) == tree.getTaxon(taxonIndex)) {
                return tip.getNumber();
            }
        }
        throw new IllegalStateException();
    }

    public int size() {
        return tips.size();
    }

    /**
     * @return the number of boolean statistics this filter does not cover (scored through the objective only)
     */
    public int getUnfilteredCount() {
        return unfiltered;
    }

    public List<Set<String>> getSets() {
        return Collections.unmodifiableList(sets);
    }

    /**
     * @return whether the move would break a constraint: clade(P) is one of them
     */
    public boolean violates(Tree tree, RootedNniMove move) {
        if (tips.isEmpty()) {
            return false;
        }
        return tips.contains(clade(tree, tree.getNode(move.getP())));
    }

    /**
     * The SPR rule (plan 7 §5.7, corrected): with L the lowest common ancestor of P and pj in the tree before the move,
     * the clades on the path from Q up to L lose the subtree, those on the path from pj up to L gain it (L excluded on
     * both), P's becomes the subtree plus j's clade after the move, and no other clade changes. The move breaks a
     * constraint iff it is one of these clades before the move and none of them after: when j is Q, for instance, Q
     * loses the subtree but P's new clade is Q's old one.
     *
     * @param clades the clades of the tree as it is, by node number (clades(tree))
     */
    public boolean violates(Tree tree, ClampedSprMove move, BitSet[] clades) {
        if (tips.isEmpty()) {
            return false;
        }
        final NodeRef p = tree.getNode(move.getP());
        final NodeRef pj = tree.getNode(move.getPj());
        final BitSet subtree = clades[move.getI()];
        final BitSet aboveP = new BitSet(tree.getNodeCount());
        for (NodeRef x = p; x != null; x = tree.getParent(x)) {
            aboveP.set(x.getNumber());
        }
        NodeRef lca = pj;
        while (!aboveP.get(lca.getNumber())) {
            lca = tree.getParent(lca);
        }

        final List<BitSet> before = new ArrayList<BitSet>();
        final List<BitSet> after = new ArrayList<BitSet>();
        BitSet jAfter = clades[move.getJ()];
        before.add(clades[p.getNumber()]);
        // (L is P itself when j is below C: then no clade above P changes)
        for (NodeRef x = tree.getParent(p); lca != p && x != lca; x = tree.getParent(x)) {
            before.add(clades[x.getNumber()]);
            BitSet lost = (BitSet) clades[x.getNumber()].clone();
            lost.andNot(subtree);
            after.add(lost);
            if (x.getNumber() == move.getJ()) {
                jAfter = lost;
            }
        }
        for (NodeRef x = pj; x != lca; x = tree.getParent(x)) {
            before.add(clades[x.getNumber()]);
            BitSet gained = (BitSet) clades[x.getNumber()].clone();
            gained.or(subtree);
            after.add(gained);
        }
        BitSet pAfter = (BitSet) jAfter.clone();
        pAfter.or(subtree);
        after.add(pAfter);

        for (BitSet clade : before) {
            if (tips.contains(clade) && !after.contains(clade)) {
                return true;
            }
        }
        return false;
    }

    /**
     * The swap rule: with L the lowest common ancestor of iP and jP, the clades on the path from iP up to L lose
     * subtree(i) and gain subtree(j), those from jP up to L the reverse (L excluded), and no other clade changes. The
     * swap breaks a constraint iff it is one of these clades before and none of them after.
     */
    public boolean violates(Tree tree, SubtreeSwapMove move, BitSet[] clades) {
        if (tips.isEmpty()) {
            return false;
        }
        final NodeRef iP = tree.getNode(move.getIP());
        final NodeRef jP = tree.getNode(move.getJP());
        final BitSet subtreeI = clades[move.getI()];
        final BitSet subtreeJ = clades[move.getJ()];
        final BitSet aboveI = new BitSet(tree.getNodeCount());
        for (NodeRef x = iP; x != null; x = tree.getParent(x)) {
            aboveI.set(x.getNumber());
        }
        NodeRef lca = jP;
        while (!aboveI.get(lca.getNumber())) {
            lca = tree.getParent(lca);
        }
        final List<BitSet> before = new ArrayList<BitSet>();
        final List<BitSet> after = new ArrayList<BitSet>();
        for (NodeRef x = iP; x != lca; x = tree.getParent(x)) {
            before.add(clades[x.getNumber()]);
            after.add(swapped(clades[x.getNumber()], subtreeI, subtreeJ));
        }
        for (NodeRef x = jP; x != lca; x = tree.getParent(x)) {
            before.add(clades[x.getNumber()]);
            after.add(swapped(clades[x.getNumber()], subtreeJ, subtreeI));
        }
        for (BitSet clade : before) {
            if (tips.contains(clade) && !after.contains(clade)) {
                return true;
            }
        }
        return false;
    }

    /**
     * The root-move rule. Above the root: the clades from Q up to the root (included) lose subtree(i), P's old clade
     * goes and P's new one is all tips. At the root: the clades from k up to C (included) gain subtree(i), and the old
     * root's clade, all tips, becomes subtree(i) plus subtree(j). The move breaks a constraint iff it is one of these
     * clades before and none of them after.
     */
    public boolean violates(Tree tree, RootChangingSprMove move, BitSet[] clades) {
        if (tips.isEmpty()) {
            return false;
        }
        final BitSet subtree = clades[move.getI()];
        final List<BitSet> before = new ArrayList<BitSet>();
        final List<BitSet> after = new ArrayList<BitSet>();
        if (move.getKind() == RootChangingSprMove.Kind.ABOVE_ROOT) {
            before.add(clades[move.getP()]);
            for (NodeRef x = tree.getNode(move.getQ()); x != null; x = tree.getParent(x)) {
                before.add(clades[x.getNumber()]);
                BitSet lost = (BitSet) clades[x.getNumber()].clone();
                lost.andNot(subtree);
                after.add(lost);
            }
        } else {
            for (NodeRef x = tree.getNode(move.getK()); x.getNumber() != move.getP(); x = tree.getParent(x)) {
                before.add(clades[x.getNumber()]);
                BitSet gained = (BitSet) clades[x.getNumber()].clone();
                gained.or(subtree);
                after.add(gained);
            }
            BitSet pAfter = (BitSet) clades[move.getJ()].clone();
            pAfter.or(subtree);
            after.add(pAfter);
        }
        for (BitSet clade : before) {
            if (tips.contains(clade) && !after.contains(clade)) {
                return true;
            }
        }
        return false;
    }

    private static BitSet swapped(BitSet clade, BitSet out, BitSet in) {
        BitSet result = (BitSet) clade.clone();
        result.andNot(out);
        result.or(in);
        return result;
    }

    /**
     * @return the tip numbers below every node, by node number
     */
    public static BitSet[] clades(Tree tree) {
        BitSet[] clades = new BitSet[tree.getNodeCount()];
        fill(tree, tree.getRoot(), clades);
        return clades;
    }

    private static BitSet fill(Tree tree, NodeRef node, BitSet[] clades) {
        BitSet bits = new BitSet(tree.getExternalNodeCount());
        if (tree.isExternal(node)) {
            bits.set(node.getNumber());
        } else {
            for (int c = 0; c < tree.getChildCount(node); ++c) {
                bits.or(fill(tree, tree.getChild(node, c), clades));
            }
        }
        clades[node.getNumber()] = bits;
        return bits;
    }

    /**
     * @return whether the node's clade is a constraint (a move that removes it breaks the constraint)
     */
    public boolean isConstrained(Tree tree, int node) {
        return !tips.isEmpty() && tips.contains(clade(tree, tree.getNode(node)));
    }

    /**
     * @return whether every constraint is a clade of the tree
     */
    public boolean areSatisfied(Tree tree) {
        for (Set<String> set : sets) {
            if (!TreeUtils.isMonophyletic(tree, set)) {
                return false;
            }
        }
        return true;
    }

    // the tip numbers below the node
    static BitSet clade(Tree tree, NodeRef node) {
        BitSet bits = new BitSet(tree.getExternalNodeCount());
        addTips(tree, node, bits);
        return bits;
    }

    private static void addTips(Tree tree, NodeRef node, BitSet bits) {
        if (tree.isExternal(node)) {
            bits.set(node.getNumber());
        } else {
            for (int c = 0; c < tree.getChildCount(node); ++c) {
                addTips(tree, tree.getChild(node, c), bits);
            }
        }
    }
}
