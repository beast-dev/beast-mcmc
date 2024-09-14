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

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

/**
 * @author Marc A Suchard
 * @author Philippe Lemey
 */
public class RecombinationHeightBoundStatistic extends TreeStatistic {

    public RecombinationHeightBoundStatistic(String name, Tree tree,
                                             List<TaxonList> taxa,
                                             boolean absoluteTime) throws TreeUtils.MissingTaxonException {

        super(name);
        this.tree = tree;

        this.leafSets = new ArrayList<>();
        for (TaxonList taxon : taxa) {
            leafSets.add(TreeUtils.getLeavesForTaxa(tree, taxon));
        }

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

        double max = Double.NEGATIVE_INFINITY;

        for (Set<String> leafSet : leafSets) {
            NodeRef recombinant = TreeUtils.getCommonAncestorNode(tree, leafSet);
            if (tree.isRoot(recombinant)) {
                return tree.getNodeHeight(tree.getRoot());
            }

            max = Math.max(max, tree.getNodeHeight(tree.getParent(recombinant)));
        }

        return max;
    }

    private double getMinHeight() {

        double min = Double.POSITIVE_INFINITY;

        for (Set<String> leafSet : leafSets) {
            NodeRef recombinant = TreeUtils.getCommonAncestorNode(tree, leafSet);

            min = Math.min(min, tree.getNodeHeight(recombinant));
        }

        return min;
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

    private final List<Set<String>> leafSets;
    private final double mostRecentTipTime;
    private final boolean isBackwards;
}
