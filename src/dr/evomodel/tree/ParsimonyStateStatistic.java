/*
 * ParsimonyStateStatistic.java
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
import dr.evolution.util.TaxonList;
import dr.inference.model.Statistic;

import java.util.Set;


/**
 * A statistic that reconstructs the parsimony state at a mrca
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: ParsimonyStateStatistic.java,v 1.11 2005/07/11 14:06:25 rambaut Exp $
 */
public class ParsimonyStateStatistic extends TreeStatistic {

    public ParsimonyStateStatistic(String name, Tree tree, TaxonList stateTaxa, TaxonList mrcaTaxa) throws TreeUtils.MissingTaxonException {

        super(name);
        this.tree = tree;
        this.stateLeafSet = TreeUtils.getLeavesForTaxa(tree, stateTaxa);
        if (mrcaTaxa != null) {
            this.mrcaLeafSet = TreeUtils.getLeavesForTaxa(tree, mrcaTaxa);
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
     * @return the parsimony tree length of the character.
     */
    public double getStatisticValue(int dim) {

        NodeRef node;

        if (mrcaLeafSet != null) {
            node = TreeUtils.getCommonAncestorNode(tree, mrcaLeafSet);
        } else {
            node = tree.getRoot();
        }
        return TreeUtils.getParsimonyState(tree, node, stateLeafSet);
    }

    private Tree tree = null;
    private Set stateLeafSet = null;
    private Set<String> mrcaLeafSet = null;
}
