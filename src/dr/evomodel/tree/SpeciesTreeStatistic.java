/*
 * SpeciesTreeStatistic.java
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
import dr.inference.model.BooleanStatistic;

import java.util.HashSet;
import java.util.Set;

/**
 * Interface for a phylogenetic tree that may contain population level data
 *
 * @author Alexei Drummond
 * @version $Id: SpeciesTreeStatistic.java,v 1.14 2005/07/11 14:06:25 rambaut Exp $
 */
public class SpeciesTreeStatistic extends TreeStatistic implements BooleanStatistic {

    public SpeciesTreeStatistic(String name, Tree speciesTree, Tree populationTree) {

        super(name);
        this.speciesTree = speciesTree;
        this.popTree = populationTree;
    }

    public void setTree(Tree tree) {
        this.popTree = tree;
    }

    public Tree getTree() {
        return popTree;
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
     * @return true if the population tree is compatible with the species tree
     */
    public boolean getBoolean(int dim) {

        if (popTree.getNodeHeight(popTree.getRoot()) < speciesTree.getNodeHeight(speciesTree.getRoot())) {
            return false;
        }

        return isCompatible(popTree.getRoot(), null);
    }


    private boolean isCompatible(NodeRef popNode, Set<String> species) {

        //int n = popNode.getNumber() - popTree.getExternalNodeCount();

        if (popTree.isExternal(popNode)) {
            Taxon speciesTaxon = (Taxon) popTree.getTaxonAttribute(popNode.getNumber(), "species");
            species.add(speciesTaxon.getId());
        } else {

            Set<String> speciesTaxa = new HashSet<String>();

            int childCount = popTree.getChildCount(popNode);
            for (int i = 0; i < childCount; i++) {
                if (!isCompatible(popTree.getChild(popNode, i), speciesTaxa)) {
                    return false;
                }
            }

            if (species != null) {
                species.addAll(speciesTaxa);

                NodeRef speciesNode = TreeUtils.getCommonAncestorNode(speciesTree, speciesTaxa);
                if (popTree.getNodeHeight(popNode) < speciesTree.getNodeHeight(speciesNode)) {
                    return false;
                }
            }
        }

        return true;
    }

    private Tree speciesTree;
    private Tree popTree;
}

