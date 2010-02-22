/*
 * ParsimonyStatistic.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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
public class ParsimonyStatistic extends Statistic.Abstract implements TreeStatistic {

    public ParsimonyStatistic(String name, Tree tree, TaxonList taxa) throws Tree.MissingTaxonException {

        super(name);
        this.tree = tree;
        this.leafSet = Tree.Utils.getLeavesForTaxa(tree, taxa);
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

        return Tree.Utils.getParsimonySteps(tree, leafSet);
    }

    private Tree tree = null;
    private Set leafSet = null;

}
