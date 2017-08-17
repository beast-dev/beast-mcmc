/*
 * ParsimonyStatistic.java
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

import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.TaxonList;
import dr.inference.model.Statistic;

import java.util.Set;


/**
 * A model component for trees
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: ParsimonyStatistic.java,v 1.13 2005/07/11 14:06:25 rambaut Exp $
 */
public class ParsimonyStatistic extends TreeStatistic {

    public ParsimonyStatistic(String name, Tree tree, TaxonList taxa) throws TreeUtils.MissingTaxonException {

        super(name);
        this.tree = tree;
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

    /**
     * @return the parsimony tree length of the character.
     */
    public double getStatisticValue(int dim) {

        return TreeUtils.getParsimonySteps(tree, leafSet);
    }

    private Tree tree = null;
    private Set leafSet = null;

}
