/*
 * TopologyHash.java
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

import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;

/**
 * A 64-bit hash of a rooted topology: FNV-1a of the canonical rooted Newick string (taxon names only, children in a
 * fixed order), so trees with the same rooted topology hash alike whatever their node numbers, child order or heights.
 *
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public final class TopologyHash {

    private TopologyHash() {
    }

    public static long of(Tree tree) {
        return fnv1a(TreeUtils.uniqueNewick(tree, tree.getRoot()));
    }

    public static String toHex(long hash) {
        return String.format("%016x", hash);
    }

    static long fnv1a(String text) {
        long hash = 0xcbf29ce484222325L;
        for (int i = 0; i < text.length(); ++i) {
            hash ^= text.charAt(i);
            hash *= 0x100000001b3L;
        }
        return hash;
    }
}
