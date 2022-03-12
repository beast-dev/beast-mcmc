/*
 * RandomEffectsTreeTraitProvider.java
 *
 * Copyright (c) 2002-2022 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.branchratemodel;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evolution.tree.TreeTraitProvider;

/**
 * @author Marc A. Suchard
 * @author Philippe Lemey
 */
public class RandomEffectsTreeTraitProvider implements TreeTraitProvider, TreeTrait<Double> {

    public RandomEffectsTreeTraitProvider(ArbitraryBranchRates branchRates,
                                          String traitName) {
        this.branchRates = branchRates;
        this.traitName = traitName;
    }

    public String getTraitName() {
        return traitName;
    }

    public Intent getIntent() {
        return Intent.BRANCH;
    }

    public TreeTrait getTreeTrait(final String key) {
        return this;
    }

    public TreeTrait[] getTreeTraits() {
        return new TreeTrait[]{this};
    }

    public Class getTraitClass() {
        return Double.class;
    }

    public boolean getLoggable() {
        return true;
    }

    public Double getTrait(final Tree tree, final NodeRef node) {
        return branchRates.getUntransformedBranchRate(tree, node);
    }

    public String getTraitString(final Tree tree, final NodeRef node) {
        return Double.toString(getTrait(tree, node));
    }

    private final ArbitraryBranchRates branchRates;
    private final String traitName;
}