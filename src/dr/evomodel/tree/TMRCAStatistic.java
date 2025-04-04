/*
 * TMRCAStatistic.java
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
 * A statistic that tracks the time of MRCA of a set of taxa
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
public class TMRCAStatistic extends TreeStatistic {

    public TMRCAStatistic(String name, Tree tree, TaxonList taxa, boolean absoluteTime, boolean forParent)
            throws TreeUtils.MissingTaxonException {
        super(name);
        this.tree = tree;
        if (absoluteTime && Taxon.getMostRecentDate() != null) {
            isBackwards = Taxon.getMostRecentDate().isBackwards();
            mostRecentTipTime = Taxon.getMostRecentDate().getAbsoluteTimeValue();
        } else {
            // give node heights or taxa don't have dates
            mostRecentTipTime = Double.NaN;
            isBackwards = false;
        }
        if (taxa != null) {
            this.leafSet = TreeUtils.getLeavesForTaxa(tree, taxa);
        } else {
            // if no taxa are given then use the root of the tree
            this.leafSet = null;
        }
        this.forParent = forParent;
    }

    public void setTree(Tree tree) {
        this.tree = tree;
    }

    public Tree getTree() {
        return tree;
    }

    public Set<String> getLeafSet() {
        return leafSet;
    }

    public int getDimension() {
        return 1;
    }

    /**
     * @return the height of the MRCA node.
     */
    public double getStatisticValue(int dim) {

        NodeRef node;
        if (leafSet != null) {
            node = TreeUtils.getCommonAncestorNode(tree, leafSet);
            if (forParent && !tree.isRoot(node)) {
                node = tree.getParent(node);
            }
        } else {
            node = tree.getRoot();
        }
        if (node == null) throw new RuntimeException("No node found that is MRCA of " + leafSet);

        if (!Double.isNaN(mostRecentTipTime)) {
            if (isBackwards) {
                return mostRecentTipTime + tree.getNodeHeight(node);
            } else {
                return mostRecentTipTime - tree.getNodeHeight(node);
            }
        } else {
            return tree.getNodeHeight(node);
        }
    }

    private Tree tree = null;
    private Set<String> leafSet = null;
    private final double mostRecentTipTime;
    private final boolean isBackwards;
    private final boolean forParent;

}
