/*
 * BranchSpecificRateBranchModel.java
 *
 * Copyright (c) 2002-2018 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evolution.tree.TreeTraitProvider;
import dr.evomodel.substmodel.BranchSpecificSubstitutionModelProvider;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.TreeParameterModel;
import dr.inference.model.*;

import java.util.List;

/**
 * @author Marc Suchard
 * @author Xiang Ji
 */
public class ArbitrarySubstitutionParameterBranchModel extends AbstractModel implements BranchModel, TreeTraitProvider{

    private final BranchSpecificSubstitutionModelProvider substitutionModelProvider;
    protected final List<CompoundParameter> substitutionParameterList;
    private int parameterDimension;
    private final TreeModel tree;
    private final TreeParameterModel parameterIndexHelper;
    private final Helper traitProvider = new Helper();


    public ArbitrarySubstitutionParameterBranchModel(String name,
                                                     BranchSpecificSubstitutionModelProvider substitutionModelProvider,
                                                     List<CompoundParameter> substitutionParameterList,
                                                     TreeModel tree) {
        super(name);
        this.substitutionModelProvider = substitutionModelProvider;
        this.tree = tree;
        this.parameterIndexHelper = new TreeParameterModel(tree,  substitutionParameterList.get(0), true);

        for (SubstitutionModel substitutionModel : substitutionModelProvider.getSubstitutionModelList()) {
            addModel(substitutionModel);
        }

        this.substitutionParameterList = substitutionParameterList;
        this.parameterDimension = substitutionParameterList.get(0).getDimension();
        for (CompoundParameter substitutionParameter : substitutionParameterList) {
            addVariable(substitutionParameter);
            traitProvider.addTrait(new SubstitutionParameterTrait(substitutionParameter.getId(), substitutionParameter, parameterIndexHelper));
            assert(substitutionParameter.getDimension() == parameterDimension);
        }
    }

    public void addSubstitutionParameter(CompoundParameter substitutionParameter) {
        if (substitutionParameter.getDimension() != parameterDimension) {
            throw new RuntimeException("Wrong dimension of the new BranchSpecificSubstitutionParameter.");
        }
        this.substitutionParameterList.add(substitutionParameter);
        addVariable(substitutionParameter);
    }

    public Parameter getSubstitutionParameterForBranch(NodeRef branch, CompoundParameter branchParameter) {
        if (!substitutionParameterList.contains(branchParameter)) {
            throw new RuntimeException("The branch parameter is not in the substitution parameter list.");
        }
        return branchParameter.getParameter(parameterIndexHelper.getParameterIndexFromNodeNumber(branch.getNumber()));
    }

    @Override
    public Mapping getBranchModelMapping(NodeRef branch) {
        return substitutionModelProvider.getBranchModelMapping(branch);
    }

    @Override
    public List<SubstitutionModel> getSubstitutionModels() {
        return substitutionModelProvider.getSubstitutionModelList();
    }

    @Override
    public SubstitutionModel getRootSubstitutionModel() {
        return substitutionModelProvider.getRootSubstitutionModel();
    }

    public SubstitutionModel getSubstitutionModelForBranch(NodeRef branch) {
        return substitutionModelProvider.getSubstitutionModel(tree, branch);
    }

    @Override
    public FrequencyModel getRootFrequencyModel() {
        return getRootSubstitutionModel().getFrequencyModel();
    }

    @Override
    public boolean requiresMatrixConvolution() {
        return false;
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        fireModelChanged(object, index);
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

    }

    @Override
    public TreeTrait[] getTreeTraits() {
        return traitProvider.getTreeTraits();
    }

    @Override
    public TreeTrait getTreeTrait(String key) {
        return traitProvider.getTreeTrait(key);
    }

    private class SubstitutionParameterTrait implements TreeTrait<Double> {

        private CompoundParameter substitutionParameter;
        private String traitName;
        private TreeParameterModel parameterIndexHelper;

        private SubstitutionParameterTrait(String name,
                                           CompoundParameter substitutionParameter,
                                           TreeParameterModel parameterIndexHelper) {
            this.substitutionParameter = substitutionParameter;
            this.traitName = name;
            this.parameterIndexHelper = parameterIndexHelper;
        }

        @Override
        public String getTraitName() {
            return traitName;
        }

        @Override
        public Intent getIntent() {
            return Intent.BRANCH;
        }

        @Override
        public Class getTraitClass() {
            return CompoundParameter.class;
        }

        @Override
        public Double getTrait(Tree tree, NodeRef node) {
            return substitutionParameter.getParameterValue(0, parameterIndexHelper.getParameterIndexFromNodeNumber(node.getNumber()));
        }

        @Override
        public String getTraitString(Tree tree, NodeRef node) {
            return getTrait(tree, node).toString();
        }

        @Override
        public boolean getLoggable() {
            return true;
        }
    }
}
