/*
 * CompoundBranchRateModel.java
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

package dr.evomodel.branchratemodel;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodelxml.branchratemodel.CompoundBranchRateModelParser;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

/**
 * Takes a collection of BranchRateModels and returns the product of the respective rates. In order for this to
 * work, one model should drive the actual rate of evolution and the others should be set up to provide
 * relative rates.
 * @author Andrew Rambaut
 * @version $Id:=$
 */
public class CompoundBranchRateModel extends AbstractBranchRateModel {

    private final List<BranchRateModel> branchRateModels = new ArrayList<BranchRateModel>();

    public CompoundBranchRateModel(Collection<BranchRateModel> branchRateModels) {
        super(CompoundBranchRateModelParser.COMPOUND_BRANCH_RATE_MODEL);
        for (BranchRateModel branchRateModel : branchRateModels) {
            addModel(branchRateModel);
            this.branchRateModels.add(branchRateModel);
        }
    }

    public void handleModelChangedEvent(Model model, Object object, int index) {
        if (index != -1) {
            fireModelChanged(null, index);
        } else {
            fireModelChanged();
        }
    }

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
    }

    protected void storeState() {
        // nothing to do
    }

    protected void restoreState() {
        // nothing to do
    }

    protected void acceptState() {
        // nothing to do
    }

    public double getBranchRate(final Tree tree, final NodeRef node) {
        double rate = 1.0;
        for (BranchRateModel branchRateModel : branchRateModels) {
            rate *= branchRateModel.getBranchRate(tree, node);
        }
        return rate;
    }

}