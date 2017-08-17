/*
 * ExternalLengthStatistic.java
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
import dr.inference.model.Statistic;

import java.util.ArrayList;
import java.util.List;

/**
 * A statistic that tracks the time of MRCA of a set of taxa
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: TMRCAStatistic.java,v 1.21 2005/07/11 14:06:25 rambaut Exp $
 */
public class ExternalLengthStatistic extends TreeStatistic {

    public ExternalLengthStatistic(String name, Tree tree, TaxonList taxa) throws TreeUtils.MissingTaxonException {
        super(name);
        this.tree = tree;
        int m = taxa.getTaxonCount();
        int n = tree.getExternalNodeCount();

        for (int i = 0; i < m; i++) {

            Taxon taxon = taxa.getTaxon(i);
            NodeRef node = null;
            boolean found = false;
            for (int j = 0; j < n; j++) {

                node = tree.getExternalNode(j);
                if (tree.getNodeTaxon(node).getId().equals(taxon.getId())) {

                    found = true;
                    break;
                }
            }

            if (!found) {
                throw new TreeUtils.MissingTaxonException(taxon);
            }

            leafSet.add(node);
        }
    }

    public void setTree(Tree tree) {
        this.tree = tree;
    }

    public Tree getTree() {
        return tree;
    }

    public int getDimension() {
        return leafSet.size();
    }

    /**
     * @return the height of the MRCA node.
     */
    public double getStatisticValue(int dim) {
        NodeRef node = leafSet.get(dim);
        return tree.getBranchLength(node);
    }

    private Tree tree = null;
    private List<NodeRef> leafSet = new ArrayList<NodeRef>();

}