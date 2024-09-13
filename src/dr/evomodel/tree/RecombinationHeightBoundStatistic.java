/*
 * RecombinationHeightBoundStatistic.java
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

package dr.evomodel.tree;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;

import java.util.Set;

/**
 * @author Marc A Suchard
 * @author Philippe Lemey
 */
public class RecombinationHeightBoundStatistic extends TreeStatistic {

    public RecombinationHeightBoundStatistic(String name, Tree tree,
                                             TaxonList taxa0,
                                             TaxonList taxa1,
                                             boolean absoluteTime) throws TreeUtils.MissingTaxonException {

        super(name);
        this.tree = tree;

        this.leafSet0 = TreeUtils.getLeavesForTaxa(tree, taxa0);
        this.leafSet1 = TreeUtils.getLeavesForTaxa(tree, taxa1);

        if (absoluteTime && Taxon.getMostRecentDate() != null) {
            isBackwards = Taxon.getMostRecentDate().isBackwards();
            mostRecentTipTime = Taxon.getMostRecentDate().getAbsoluteTimeValue();
        } else {
            // give node heights or taxa don't have dates
            mostRecentTipTime = Double.NaN;
            isBackwards = false;
        }
    }

    public void setTree(Tree tree) {
        this.tree = tree;
    }

    public Tree getTree() {
        return tree;
    }

    public int getDimension() {
        return 2;
    }

    public double getStatisticValue(int dim) {

        final double height;
        if (dim == 0) {
            height = getMaxHeight();
        } else if (dim == 1) {
            height = getMinHeight();
        } else {
            throw new IllegalArgumentException("Unknown dimension");
        }

        return transformHeight(height);
    }

    private double getMaxHeight() {
        NodeRef recombinant0 = TreeUtils.getCommonAncestorNode(tree, leafSet0);
        NodeRef recombinant1 = TreeUtils.getCommonAncestorNode(tree, leafSet1);

        if (tree.isRoot(recombinant0) || tree.isRoot(recombinant1)) {
            return tree.getNodeHeight(tree.getRoot());
        }

        double height0 = tree.getNodeHeight(tree.getParent(recombinant0));
        double height1 = tree.getNodeHeight(tree.getParent(recombinant1));

        return Math.max(height0, height1);
    }

    private double getMinHeight() {
        NodeRef recombinant0 = TreeUtils.getCommonAncestorNode(tree, leafSet0);
        NodeRef recombinant1 = TreeUtils.getCommonAncestorNode(tree, leafSet1);

        double height0 = tree.getNodeHeight(recombinant0);
        double height1 = tree.getNodeHeight(recombinant1);

        return Math.min(height0, height1);
    }

    private double transformHeight(double height) {
        if (!Double.isNaN(mostRecentTipTime)) {
            if (isBackwards) {
                return mostRecentTipTime + height;
            } else {
                return mostRecentTipTime - height;
            }
        } else {
            return height;
        }
    }

    private Tree tree;

    private final Set<String> leafSet0;
    private final Set<String> leafSet1;

    private final double mostRecentTipTime;
    private final boolean isBackwards;
}
