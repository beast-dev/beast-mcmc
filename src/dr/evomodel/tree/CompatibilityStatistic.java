/*
 * CompatibilityStatistic.java
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
import dr.inference.model.BooleanStatistic;

import java.util.BitSet;
import java.util.HashSet;
import java.util.Set;

/**
 * Tests whether 2 (possibly unresolved) trees are compatible
 * *
 *
 * @author Andrew Rambaut
 */
public class CompatibilityStatistic extends TreeStatistic implements BooleanStatistic {

    public CompatibilityStatistic(String name, Tree tree1, Tree tree2) throws TreeUtils.MissingTaxonException {

        super(name);
        this.tree = tree1;

        intersection = new BitSet(tree1.getExternalNodeCount());
        clades = new HashSet<BitSet>();
        getClades(tree1, tree2, tree2.getRoot(), null, clades);

        for (int i = 0; i < tree1.getTaxonCount(); i++) {
            String id = tree1.getTaxonId(i);
            if (tree2.getTaxonIndex(id) == -1) {
                throw new TreeUtils.MissingTaxonException(tree1.getTaxon(i));
            }
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
     * @return boolean result of test.
     */
    public double getStatisticValue(int dim) {
        return getBoolean(dim) ? 1.0 : 0.0;
    }

    /**
     * @return boolean result of test.
     */
    public boolean getBoolean(int dim) {
        return isCompatible(tree, tree.getRoot(), null);
    }

    private boolean isCompatible(Tree tree, NodeRef node, BitSet leaves) {
        if (tree.isExternal(node)) {
            leaves.set(node.getNumber());
            return true;
        } else {

            BitSet ls = new BitSet(tree.getExternalNodeCount());

            for (int i = 0; i < tree.getChildCount(node); i++) {

                NodeRef node1 = tree.getChild(node, i);

                if (!isCompatible(tree, node1, ls)) {
                    // as soon as we have an incompatibility break out...
                    return false;
                }
            }

            if (leaves != null) {
                if (!clades.contains(ls)) {
                    // except for the root clade...
                    for (BitSet clade : clades) {
                        intersection.clear();
                        intersection.or(clade);
                        intersection.and(ls);

                        int card = intersection.cardinality();
                        if (card != 0 &&
                                card != ls.cardinality() &&
                                card != clade.cardinality()) {
                            return false;
                        }
                    }
                }

                leaves.or(ls);
            }
        }
        return true;
    }

    private void getClades(Tree referenceTree, Tree tree, NodeRef node, BitSet leaves, Set<BitSet> clades) {

        if (tree.isExternal(node)) {
            String taxonId = tree.getNodeTaxon(node).getId();
            for (int i = 0; i < referenceTree.getExternalNodeCount(); i++) {
                NodeRef n = referenceTree.getExternalNode(i);

                if (taxonId.equals(referenceTree.getNodeTaxon(n).getId())) {
                    leaves.set(n.getNumber());
                }
            }
        } else {

            BitSet ls = new BitSet(tree.getExternalNodeCount());

            for (int i = 0; i < tree.getChildCount(node); i++) {

                NodeRef node1 = tree.getChild(node, i);

                getClades(referenceTree, tree, node1, ls, clades);
            }

            if (leaves != null) {
                // except for the root clade...
                leaves.or(ls);
                clades.add(ls);
            }

        }
    }

    private Tree tree;
    private final Set<BitSet> clades;
    private final BitSet intersection;

}