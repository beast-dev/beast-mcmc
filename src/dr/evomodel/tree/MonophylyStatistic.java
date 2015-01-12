/*
 * MonophylyStatistic.java
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
import dr.inference.model.BooleanStatistic;

import java.util.Collections;
import java.util.Set;

/**
 * Performs monophyly test given a taxonList
 *
 * @author Roald Forsberg
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: MonophylyStatistic.java,v 1.16 2005/07/11 14:06:25 rambaut Exp $
 */
public class MonophylyStatistic extends BooleanStatistic implements TreeStatistic {

    public MonophylyStatistic(String name, Tree tree, TaxonList taxa, TaxonList ignore) throws Tree.MissingTaxonException {

        this(name, tree, taxa, ignore, false);

    }

    public MonophylyStatistic(String name, Tree tree, TaxonList taxa, TaxonList ignore, boolean inverse) throws Tree.MissingTaxonException {

        super(name);
        this.tree = tree;
        this.leafSet = Tree.Utils.getLeavesForTaxa(tree, taxa);
        if (ignore != null) {
            this.ignoreLeafSet = Tree.Utils.getLeavesForTaxa(tree, ignore);
        } else {
            this.ignoreLeafSet = Collections.emptySet();
        }
        this.inverse = inverse;

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
    public boolean getBoolean(int dim) {
        boolean monophyletic = Tree.Utils.isMonophyletic(this.tree, this.leafSet, this.ignoreLeafSet);
        if (inverse){
            return !monophyletic;
        } else {
            return monophyletic;
        }
    }

    private Tree tree = null;
    private Set<String> leafSet = null;
    private Set<String> ignoreLeafSet = null;
    private boolean inverse = false;

}
