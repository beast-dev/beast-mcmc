/*
 * HomogeneousBranchModel.java
 *
 * Copyright (C) 2002-2012 Alexei Drummond, Andrew Rambaut & Marc A. Suchard
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

package dr.app.beagle.evomodel.branchmodel;

import dr.app.beagle.evomodel.substmodel.SubstitutionModel;
import dr.evolution.tree.NodeRef;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Variable;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Andrew Rambaut
 * @author Filip Bielejec
 * @version $Id$
 */
public class HomogeneousBranchModel extends AbstractModel implements BranchModel{
    private final SubstitutionModel substitutionModel;

    public HomogeneousBranchModel(SubstitutionModel substitutionModel) {
        super("HomogeneousBranchModel");
        this.substitutionModel = substitutionModel;
        addModel(substitutionModel);
    }

    public Mapping getBranchModelMapping(NodeRef node) {
        return DEFAULT;
    }

    @Override
    public List<SubstitutionModel> getSubstitutionModels() {
        List<SubstitutionModel> substitutionModels = new ArrayList<SubstitutionModel>();
        substitutionModels.add(substitutionModel);
        return substitutionModels;
    }

    @Override
    public SubstitutionModel getRootSubstitutionModel() {
        return substitutionModel;
    }

    @Override
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
}
