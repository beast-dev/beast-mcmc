/*
 * HierarchicalRateAggregatedPattern.java
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

package dr.evolution.alignment;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.TreeSet;

/**
 * The patterns of a 4-state pattern list, ordered so that the patterns sharing an aggregated (reduced) character
 * tree are contiguous, for BEAGLE's blocked aggregation of the hierarchical representation.
 * <p>
 * A pattern's block is the set of internal nodes of the character tree with an observed state (a state code below
 * the state count) below them. Every rooted 4-state tree gives two reduced blocks plus the full block. The order is
 * [first reduced block][full block and all-missing patterns][second reduced block], stable within each block; the
 * first reduced block is the one whose kept-node set has the smaller bit pattern. Weights move with their patterns.
 * The list carries no site-to-pattern map, like the {@link Patterns} copies the tree data likelihood already uses.
 *
 * @author Jeff Thorne
 * @author Tae-Kun Seo
 * @author Sangchul Choi
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public class HierarchicalRateAggregatedPattern extends Patterns {

    private final int[] blockSizes = new int[3];

    /**
     * @param source  a list of 4-state patterns, without uncertain states
     * @param triples the character tree as BEAGLE triples (parent, child, child) in post-order: states are tips 0..3,
     *                internal nodes 4..6, root last (HierarchicalSubstitutionModel.getBeagleModelTriples())
     */
    public HierarchicalRateAggregatedPattern(PatternList source, double[] triples) {
        super(source.getDataType(), source);
        final int stateCount = source.getStateCount();
        if (stateCount != 4 || triples.length < 3 * (stateCount - 1) || source.areUncertain()) {
            throw new IllegalArgumentException("blocked aggregation needs 4-state patterns without uncertain states"
                    + " and a 4-state character tree");
        }
        final int nodeCount = 2 * stateCount - 1;
        final int[] parent = new int[nodeCount];
        Arrays.fill(parent, -1);
        for (int t = 0; t < stateCount - 1; t++) {
            final int m = (int) triples[3 * t];
            parent[(int) triples[3 * t + 1]] = m;
            parent[(int) triples[3 * t + 2]] = m;
        }
        final int[] keptOfState = new int[stateCount];
        for (int x = 0; x < stateCount; x++) {
            for (int m = parent[x]; m >= 0; m = parent[m]) {
                keptOfState[x] |= 1 << (m - stateCount);
            }
        }
        final int all = (1 << (stateCount - 1)) - 1;
        final int count = source.getPatternCount();
        final int[] kept = new int[count];
        final TreeSet<Integer> reduced = new TreeSet<Integer>();
        for (int k = 0; k < count; k++) {
            final int[] pattern = source.getPattern(k);
            int bits = 0;
            for (int state : pattern) {
                if (state >= 0 && state < stateCount) {
                    bits |= keptOfState[state];
                }
            }
            kept[k] = bits;
            if (bits != 0 && bits != all) {
                reduced.add(bits);
            }
        }
        final List<Integer> keys = new ArrayList<Integer>(reduced);
        if (keys.size() > 2) {
            throw new IllegalStateException("a 4-state character tree has at most two reduced blocks");
        }
        final List<Integer> order = new ArrayList<Integer>(count);
        for (int group = 0; group < 3; group++) {
            for (int k = 0; k < count; k++) {
                if (groupOf(kept[k], keys) == group) {
                    order.add(k);
                }
            }
        }
        patterns = new int[count][];
        weights = new double[count];
        for (int i = 0; i < count; i++) {
            final int k = order.get(i);
            patterns[i] = source.getPattern(k);
            weights[i] = source.getPatternWeight(k);
            blockSizes[groupOf(kept[k], keys)]++;
        }
        patternCount = count;
    }

    private static int groupOf(int bits, List<Integer> keys) {
        if (keys.size() > 0 && bits == keys.get(0)) {
            return 0;
        }
        if (keys.size() > 1 && bits == keys.get(1)) {
            return 2;
        }
        return 1;
    }

    /**
     * @return the number of patterns in the first reduced block, the full block (with all-missing patterns) and the
     * second reduced block, in list order
     */
    public int[] getBlockSizes() {
        return blockSizes.clone();
    }
}
