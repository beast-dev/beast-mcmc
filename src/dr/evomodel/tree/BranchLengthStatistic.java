/*
 * BranchLengthStatistic.java
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
 * A statistic that extracts the length of the stem branch of a set of taxa
 *
 * @author Andrew Rambaut
 */
public class BranchLengthStatistic extends TreeStatistic {

    public BranchLengthStatistic(String name, Tree tree, TaxonList taxa)
            throws TreeUtils.MissingTaxonException {
        super(name);
        this.tree = tree;
        if (taxa != null) {
            this.leafSet = TreeUtils.getLeavesForTaxa(tree, taxa);
        } else {
            throw new IllegalArgumentException("taxa cannot be null");
        }
    }

    public void setTree(Tree tree) {
        this.tree = tree;
    }

    public Tree getTree() {
        return tree;
    }

    public int getDimension() {
        return 1;
    }

    /**
     * @return the height of the MRCA node.
     */
    public double getStatisticValue(int dim) {
        NodeRef node = TreeUtils.getCommonAncestorNode(tree, leafSet);
        
        if (node == null) {
            throw new RuntimeException("No node found that is MRCA of " + leafSet);
        }

        return tree.getBranchLength(node);
    }

    private Tree tree = null;
    private Set<String> leafSet = null;

}
