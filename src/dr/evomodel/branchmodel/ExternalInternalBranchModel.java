/*
 * ExternalInternalBranchModel.java
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

package dr.evomodel.branchmodel;

import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evolution.tree.NodeRef;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Variable;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Andrew Rambaut
 * @author Filip Bielejec
 * @author Marc A. Suchard
 * @version $Id$
 */
public class ExternalInternalBranchModel extends AbstractModel implements BranchModel {
    public ExternalInternalBranchModel(TreeModel tree, SubstitutionModel externalSubstModel, SubstitutionModel internalSubstModel) {
        super("ExternalInternalBranchModel");

        this.tree = tree;
        this.externalSubstModel = externalSubstModel;
        this.internalSubstModel = internalSubstModel;

        addModel(tree);
        addModel(externalSubstModel);
        addModel(internalSubstModel);
    }

    public Mapping getBranchModelMapping(final NodeRef node) {
        return new Mapping() {
            public int[] getOrder() {
                return new int[] { tree.isExternal(node) ? 0 : 1 };
            }

            public double[] getWeights() {
                return new double[] { 1.0 };
            }
        };
    }

//    @Override // use java 1.5
    public List<SubstitutionModel> getSubstitutionModels() {
        List<SubstitutionModel> substitutionModels = new ArrayList<SubstitutionModel>();
        substitutionModels.add(externalSubstModel);
        substitutionModels.add(internalSubstModel);
        return substitutionModels;
    }

//    @Override
    public SubstitutionModel getRootSubstitutionModel() {
        return internalSubstModel;
    }

    public FrequencyModel getRootFrequencyModel() {
        return getRootSubstitutionModel().getFrequencyModel();
    }

//    @Override
    public boolean requiresMatrixConvolution() {
        return false;
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        fireModelChanged();
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
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

    private final TreeModel tree;
    private final SubstitutionModel externalSubstModel;
    private final SubstitutionModel internalSubstModel;
}
