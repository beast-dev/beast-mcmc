/*
 * CollectionAction.java
 *
 * Copyright © 2002-2026, the BEAST Development Team.
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
 */

package dr.app.tools.treeannotator;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;

import java.util.Collection;
import java.util.LinkedHashSet;
import java.util.Set;

class CollectHeightsAction implements CladeAction {
    public static final boolean STORE_TIP_HEIGHTS = true;

    public CollectHeightsAction() {
    }

    @Override
    public void actOnClade(Clade clade, Tree tree, NodeRef node) {
        collectHeightsForClade(clade, tree, node);
    }

    @Override
    public boolean expectAllClades() {
        return true;
    }

    private void collectHeightsForClade(Clade clade, Tree tree, NodeRef node) {

        if (tree.isExternal(node)) {
            assert clade.getSize() == 1;

            if (STORE_TIP_HEIGHTS) {
                clade.addHeightValue(tree.getNodeHeight(node));
            }
        } else {
            assert tree.getChildCount(node) == 2;

            clade.addHeightValue(tree.getNodeHeight(node));
        }
    }


}
