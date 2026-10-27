/*
 * TreeEdits.java
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
 * Topology edits for the search, without their own edit transaction: the caller holds one
 * ({@code beginTreeEdit}/{@code endTreeEdit}), so that an edit and the height changes that make it feasible are one
 * transaction. The operators' versions ({@code AbstractTreeOperator.exchangeNodes}) open their own, and a nested
 * {@code beginTreeEdit} throws.
 *
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public final class TreeEdits {

    private TreeEdits() {
    }

    /**
     * Swap two nodes between their parents: i moves under j's parent and j under i's parent. Neither may be the
     * other's ancestor, and their parents must differ.
     */
    public static void exchange(TreeModel tree, NodeRef i, NodeRef j) {
        final NodeRef iParent = tree.getParent(i);
        final NodeRef jParent = tree.getParent(j);
        if (iParent == null || jParent == null || iParent == jParent) {
            throw new IllegalArgumentException("exchange needs two non-root nodes with different parents");
        }
        // removeChild does not push an event; addChild flags the parent (and with it its children)
        tree.removeChild(iParent, i);
        tree.removeChild(jParent, j);
        tree.addChild(jParent, i);
        tree.addChild(iParent, j);
    }
}
