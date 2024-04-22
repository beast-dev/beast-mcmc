/*
 * MonophylyStatistic.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.tree;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.inference.model.BooleanStatistic;

import java.util.Collections;
import java.util.Set;

/**
 * Performs ancestor on stem test
 *
 * @author Marc A. Suchard
 */
public class AncestorOnStemStatistic extends TreeStatistic implements BooleanStatistic {

    public AncestorOnStemStatistic(String name,
                                   Tree tree,
                                   Taxon ancestor,
                                   TaxonList taxa) throws TreeUtils.MissingTaxonException {

        super(name);
        this.tree = tree;
        this.ancestor = ancestor;
        this.leafSet = TreeUtils.getLeavesForTaxa(tree, taxa);
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

    public double getStatisticValue(int dim) {
        return getBoolean(dim) ? 1.0 : 0.0;
    }

    public boolean getBoolean(int dim) {

        NodeRef mrca = TreeUtils.getCommonAncestorNode(tree, leafSet);
        if (mrca == tree.getRoot()) {
            return false;
        }

        NodeRef parent = tree.getParent(mrca);
        NodeRef child = getOtherChild(parent, mrca);
        Taxon taxon = tree.getNodeTaxon(child);
        if (taxon == null || taxon != ancestor) {
            return false;
        }

        double parentHeight = tree.getNodeHeight(parent);
        double grandParentHeight = tree.isRoot(parent) ?
                Double.POSITIVE_INFINITY : tree.getNodeHeight(tree.getParent(parent));

        return (grandParentHeight > parentHeight);
    }

    private NodeRef getOtherChild(NodeRef parent, NodeRef child) {
        NodeRef otherChild = tree.getChild(parent, 0);
        if (otherChild == child) {
            otherChild = tree.getChild(parent, 1);
        }
        return otherChild;
    }

    private Tree tree;
    private final Taxon ancestor;
    private final Set<String> leafSet;
}
