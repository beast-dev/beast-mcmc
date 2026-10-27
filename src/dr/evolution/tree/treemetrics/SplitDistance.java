/*
 * SplitDistance.java
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

package dr.evolution.tree.treemetrics;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;

import java.util.BitSet;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * Symmetric Robinson-Foulds distances between two trees on the same taxa, rooted (clades) and unrooted (splits), as
 * counts and normalized by their maxima for binary trees. RobinsonFouldsMetric counts the clades of the first tree
 * missing from the second only, and over all clades, tips and root included; here both directions count and only
 * the informative clades and splits do.
 * <ul>
 * <li>Rooted: the clades of the internal nodes other than the root (n - 2 in a binary tree), so the maximum is
 * 2 (n - 2).</li>
 * <li>Unrooted: the bipartitions of the internal edges, ignoring the root (n - 3 in a binary tree), so the maximum is
 * 2 (n - 3). Trees that differ only in the root position are 0 apart.</li>
 * </ul>
 * Taxa are matched by id.
 *
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public final class SplitDistance {

    private SplitDistance() {
    }

    /**
     * @return the number of informative clades in one tree and not the other, both ways
     */
    public static int rooted(Tree tree1, Tree tree2) {
        Map<String, Integer> index = taxonIndex(tree1, tree2);
        return symmetricDifference(clades(tree1, index), clades(tree2, index));
    }

    /**
     * @return the number of informative splits in one tree and not the other, both ways
     */
    public static int unrooted(Tree tree1, Tree tree2) {
        Map<String, Integer> index = taxonIndex(tree1, tree2);
        return symmetricDifference(splits(tree1, index), splits(tree2, index));
    }

    /**
     * @return rooted(tree1, tree2) / (2 (n - 2)), in [0, 1] for binary trees (0 for n < 3)
     */
    public static double normalizedRooted(Tree tree1, Tree tree2) {
        final int max = 2 * (tree1.getExternalNodeCount() - 2);
        return max > 0 ? rooted(tree1, tree2) / (double) max : 0.0;
    }

    /**
     * @return unrooted(tree1, tree2) / (2 (n - 3)), in [0, 1] for binary trees (0 for n < 4)
     */
    public static double normalizedUnrooted(Tree tree1, Tree tree2) {
        final int max = 2 * (tree1.getExternalNodeCount() - 3);
        return max > 0 ? unrooted(tree1, tree2) / (double) max : 0.0;
    }

    /**
     * @return the tip sets of the internal nodes other than the root, by the taxon index
     */
    static Set<BitSet> clades(Tree tree, Map<String, Integer> index) {
        Set<BitSet> clades = new HashSet<BitSet>();
        collect(tree, tree.getRoot(), index, clades, null);
        return clades;
    }

    /**
     * @return the bipartitions of the internal edges, each as the side without taxon 0
     */
    static Set<BitSet> splits(Tree tree, Map<String, Integer> index) {
        Set<BitSet> splits = new HashSet<BitSet>();
        collect(tree, tree.getRoot(), index, null, splits);
        return splits;
    }

    // post-order: the tips below each node; internal non-root nodes give a clade, and their edges a split unless
    // the other side is a single tip (an edge of the root's child above a tip)
    private static BitSet collect(Tree tree, NodeRef node, Map<String, Integer> index, Set<BitSet> clades,
                                  Set<BitSet> splits) {
        final BitSet below = new BitSet();
        if (tree.isExternal(node)) {
            below.set(index.get(tree.getNodeTaxon(node).getId()));
            return below;
        }
        for (int c = 0; c < tree.getChildCount(node); ++c) {
            below.or(collect(tree, tree.getChild(node, c), index, clades, splits));
        }
        if (!tree.isRoot(node)) {
            if (clades != null) {
                clades.add(below);
            }
            if (splits != null) {
                final int n = index.size();
                final int size = below.cardinality();
                if (size > 1 && size < n - 1) {
                    BitSet split = (BitSet) below.clone();
                    if (split.get(0)) {
                        split.flip(0, n);
                    }
                    splits.add(split);
                }
            }
        }
        return below;
    }

    static Map<String, Integer> taxonIndex(Tree tree1, Tree tree2) {
        Map<String, Integer> index = new HashMap<String, Integer>();
        for (int i = 0; i < tree1.getTaxonCount(); ++i) {
            index.put(tree1.getTaxonId(i), i);
        }
        boolean same = tree2.getTaxonCount() == index.size();
        for (int i = 0; same && i < tree2.getTaxonCount(); ++i) {
            same = index.containsKey(tree2.getTaxonId(i));
        }
        if (!same) {
            throw new IllegalArgumentException("The two trees have different taxa");
        }
        return index;
    }

    private static int symmetricDifference(Set<BitSet> a, Set<BitSet> b) {
        int count = 0;
        for (BitSet x : a) {
            if (!b.contains(x)) {
                ++count;
            }
        }
        for (BitSet x : b) {
            if (!a.contains(x)) {
                ++count;
            }
        }
        return count;
    }
}
