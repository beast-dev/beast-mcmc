/*
 * SetHeightsAction.java
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

import dr.evolution.tree.MutableTree;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.math.distributions.KernelDensityEstimatorDistribution;
import dr.math.distributions.NormalKDEDistribution;
import dr.util.HeapSort;

import java.util.List;

public class SetCladeHeightsAction implements CladeAction {

    // the heights to set for the root of the tree
    private final List<Double> rootHeights;

    SetCladeHeightsAction(List<Double> rootHeights) {
        this.rootHeights = rootHeights;
    }

    @Override
    public void actOnClade(Clade clade, Tree tree, NodeRef node) {
        BiClade biclade = (BiClade)clade;

        if (tree.isRoot(node)) {
            setCladeHeights(biclade, rootHeights);
            return;
        }
        if (clade.getSize() == 1) {
            return;
        }

        List<Double> heights = biclade.getHeightValues();

        if (clade.getSize() > 1) {
            List<Double> leftHeights = biclade.getHeightValues();
            List<Double> rightHeights = biclade.getHeightValues();
        }


    }

    @Override
    public boolean expectAllClades() {
        return true;
    }

    public void setCladeHeights(BiClade clade, List<Double> heights) {
    }
}
