/*
 * RLTVLoggerOnTree.java
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

package dr.evomodel.tree.randomlocalmodel;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;

/**
 * @author Marc A. Suchard
 */
public class RLTVLoggerOnTree implements TreeTrait<Double> {

    public static final String TRAIT_NAME = "changed";

    public RLTVLoggerOnTree(RandomLocalTreeVariable treeVariable) {
        this.treeVariable = treeVariable;
    }

    public String getTraitName() {
        return TRAIT_NAME;
    }

    public Intent getIntent() {
        return Intent.BRANCH;
    }

    public Class getTraitClass() {
        return Double.class;
    }

    public boolean getLoggable() {
        return true;
    }

    public Double getTrait(final Tree tree, final NodeRef node) {
        return treeVariable.isVariableSelected(tree, node) ? 1.0 : 0.0;
    }

    public String getTraitString(final Tree tree, final NodeRef node) {
        return Double.toString(getTrait(tree, node));
    }

    private RandomLocalTreeVariable treeVariable;
}
