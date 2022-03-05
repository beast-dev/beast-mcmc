/*
 * BranchRateModel.java
 *
 * Copyright (c) 2002-2020 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.evolution.tree.*;
import dr.evomodel.tree.AncestralTraitTreeModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

/**
 * @author Marc A. Suchard
 */
public class AncestralTraitBranchRateModel extends AbstractBranchRateModel {

    private final BranchRateModel branchRates;
    private final AncestralTraitTreeModel treeModel;
    private final Tree originalTree;

    public AncestralTraitBranchRateModel(BranchRateModel branchRates, AncestralTraitTreeModel treeModel) {
        super( "ancestral." + branchRates.getModelName());
        this.branchRates = branchRates;
        this.treeModel = treeModel;
        this.originalTree = treeModel.getOriginalTree();

        addModel(branchRates);
    }

    @Override
    public double getBranchRate(Tree tree, NodeRef node) {

        assert (tree == treeModel);

        NodeRef originalNode = treeModel.getOriginalNode(node);
        if (originalNode != null) {
            return branchRates.getBranchRate(originalTree, originalNode);
        } else {
            return 1.0;
        }
    }
    
    @Override
    public Mapping getBranchRateModelMapping(final Tree tree, final NodeRef node) {
        
        return new Mapping() {
            public double[] getRates() {
                return new double[] { getBranchRate(tree, node) };
            }

            public double[] getWeights() {
                return new double[] { 1.0 };
            }
        };
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == branchRates) {
            fireModelChanged(model);
        } else {
            throw new IllegalArgumentException("Unknown model");
        }
    }

    @Override
    protected void storeState() {

    }

    @Override
    protected void restoreState() {

    }

    @Override
    protected void acceptState() {

    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        throw new IllegalArgumentException("Unknown variable");
    }
}
