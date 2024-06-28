/*
 * MultiplicativeBranchRateModel.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.evomodel.branchratemodel;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodelxml.branchratemodel.MultiplicativeBranchRateModelParser;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

import java.util.List;

/**
 * @author Xiang Ji
 * @author Karthik Gangavarapu
 * @author Marc Suchard
 */
public class MultiplicativeBranchRateModel extends AbstractBranchRateModel {

    private final List<AbstractBranchRateModel> branchRateModels;
    private final int dim;

    public MultiplicativeBranchRateModel(List<AbstractBranchRateModel> branchRateModels) {
        super(MultiplicativeBranchRateModelParser.MULTIPLICATIVE_BRANCH_RATES);

        this.branchRateModels = branchRateModels;
        this.dim = branchRateModels.size();

        for (AbstractBranchRateModel branchRateModel : branchRateModels) {
            addModel(branchRateModel);
        }

        if (dim != 2) {
            throw new RuntimeException("Not yet tested.  Be careful with underflow errors.");
        }
    }

    @Override
    public double getBranchRate(Tree tree, NodeRef node) {
        double rate = branchRateModels.get(0).getBranchRate(tree, node);
        for (int i = 1; i < dim; i++) {
            rate *= branchRateModels.get(i).getBranchRate(tree, node);
        }
        return rate;
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        fireModelChanged();
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
        fireModelChanged();
    }
}
