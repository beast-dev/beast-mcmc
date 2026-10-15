/*
 * SubtreeMoveHook.java
 *
 * Copyright © 2026 the BEAST Development Team
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

package dr.evomodel.operators;

import dr.evolution.tree.NodeRef;

/**
 * Lets a subtree prune-and-regraft operator propose changes to state attached to the tree
 * (e.g. ancestral states) jointly with its topology change.
 *
 * The subtree below {@code node} is pruned from beside {@code sibling} (its parent node goes with it)
 * and regrafted on the branch above {@code destination}.
 */
public interface SubtreeMoveHook {

    /**
     * Called before the tree is edited, while it still has its original topology and heights.
     */
    void beforeSubtreeMove(NodeRef node, NodeRef sibling, NodeRef destination);

    /**
     * Called after the tree has been edited.
     *
     * @return the log Hastings ratio contributed by the hook's own proposal
     */
    double afterSubtreeMove();
}
