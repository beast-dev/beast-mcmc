/*
 * TimeTreeMove.java
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

import dr.evomodel.tree.TreeModel;

import java.util.LinkedHashSet;
import java.util.Set;

/**
 * Two moves applied one after the other, the second built for the tree the first leaves (plan 7 §5.2: two root-child
 * NNIs move the root two edges). The local blocks are the union of both moves' blocks after the move.
 *
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public final class ChainedMove implements TimeTreeMove {

    private final TimeTreeMove first;
    private final TimeTreeMove second;

    public ChainedMove(TimeTreeMove first, TimeTreeMove second) {
        this.first = first;
        this.second = second;
    }

    public TimeTreeMove getFirst() {
        return first;
    }

    public TimeTreeMove getSecond() {
        return second;
    }

    // the second move can only be checked once the first is applied
    @Override
    public boolean isValid(TreeModel tree) {
        return first.isValid(tree);
    }

    @Override
    public void apply(TreeModel tree, BranchRateWarmStart warmStart) {
        first.apply(tree, warmStart);
        second.apply(tree, warmStart);
    }

    @Override
    public void applyTopology(TreeModel tree) {
        first.applyTopology(tree);
        second.applyTopology(tree);
    }

    @Override
    public boolean needsRepair(TreeModel tree) {
        return first.needsRepair(tree);
    }

    @Override
    public int[] localHeights(TreeModel tree) {
        return union(first.localHeights(tree), second.localHeights(tree));
    }

    @Override
    public int[] localRates(TreeModel tree) {
        return union(first.localRates(tree), second.localRates(tree));
    }

    @Override
    public int[] touchedNodes() {
        return union(first.touchedNodes(), second.touchedNodes());
    }

    @Override
    public String describe() {
        return first.describe() + " + " + second.describe();
    }

    @Override
    public String toString() {
        return describe();
    }

    private static int[] union(int[] x, int[] y) {
        Set<Integer> set = new LinkedHashSet<Integer>();
        for (int i : x) {
            set.add(i);
        }
        for (int i : y) {
            set.add(i);
        }
        int[] result = new int[set.size()];
        int k = 0;
        for (int i : set) {
            result[k++] = i;
        }
        return result;
    }
}
