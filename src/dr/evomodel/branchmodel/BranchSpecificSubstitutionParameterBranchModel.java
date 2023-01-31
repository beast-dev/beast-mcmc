/*
 * BranchSpecificRateBranchModel.java
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

package dr.evomodel.branchmodel;

import dr.evolution.tree.NodeRef;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.DifferentiableBranchRates;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.ParameterReplaceableSubstitutionModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.*;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author Marc Suchard
 * @author Xiang Ji
 */
public class BranchSpecificSubstitutionParameterBranchModel extends AbstractModel implements BranchModel {

    private final ParameterReplaceableSubstitutionModel substitutionModel;
    private final TreeModel tree;
    private final List<SubstitutionModel> substitutionModelList;
    private Map<BranchRateModel, CompoundParameter> substitutionParameterMap = new HashMap<>();

    public BranchSpecificSubstitutionParameterBranchModel(String name,
                                                          List<Parameter> substitutionParameterList,
                                                          List<BranchRateModel> branchRateModelList,
                                                          ParameterReplaceableSubstitutionModel substitutionModel,
                                                          TreeModel tree) {
        super(name);
        this.substitutionModel = substitutionModel;
        this.tree = tree;
        if (substitutionParameterList.size() != branchRateModelList.size()) {
            throw new RuntimeException("Dimension mismatch!");
        }
        this.substitutionModelList = constructSubstitutionModels(substitutionParameterList, branchRateModelList);
        for (BranchRateModel branchRateModel : branchRateModelList) {
            addModel(branchRateModel);
        }
    }

    private List<SubstitutionModel> constructSubstitutionModels(List<Parameter> substitutionParameterList,
                                                                List<BranchRateModel> branchRateModelList) {
        List<SubstitutionModel> substitutionModelList = new ArrayList<>();
        for (int i = 0; i < tree.getNodeCount(); i++) {
            List<Parameter> newSubstitutionParameterList = new ArrayList<>();
            for (int j = 0; j < branchRateModelList.size(); j++) {
                Parameter rootParameter = substitutionParameterList.get(j);
                Parameter branchProxyParameter = new ProxySubstitutionParameter(branchRateModelList.get(j), tree, rootParameter, i);
                mapParamter(branchRateModelList.get(j), branchProxyParameter);
                newSubstitutionParameterList.add(branchProxyParameter);
            }
            substitutionModelList.add(substitutionModel.factory(substitutionParameterList, newSubstitutionParameterList));
        }
        return substitutionModelList;
    }

    public CompoundParameter getBranchSpecificParameters(DifferentiableBranchRates branchRateModel) {
        return substitutionParameterMap.get(branchRateModel);
    }

    private void mapParamter(BranchRateModel branchRateModel, Parameter branchProxyParameter) {
        CompoundParameter parameters = substitutionParameterMap.get(branchRateModel);
        if (parameters == null) {
            parameters = new CompoundParameter("branchSpecific." + branchRateModel.getId());
            substitutionParameterMap.put(branchRateModel, parameters);
        }
        parameters.addParameter(branchProxyParameter);
    }

    private class ProxySubstitutionParameter extends Parameter.Proxy implements ModelListener {

        private final BranchRateModel branchRateModel;
        private final int nodeNum;
        private final TreeModel tree;
        private final Parameter rootParameter;

        ProxySubstitutionParameter(BranchRateModel branchRateModel,
                                   TreeModel tree,
                                   Parameter rootParameter,
                                   int nodeNum) {
            super(null, 1);
            this.branchRateModel = branchRateModel;
            this.nodeNum = nodeNum;
            this.tree = tree;
            this.rootParameter = rootParameter;
            branchRateModel.addModelListener(this);
            addVariable(rootParameter);
        }

        @Override
        public double getParameterValue(int dim) {
            assert(dim == 0);
            NodeRef node = tree.getNode(nodeNum);
            if (tree.isRoot(node)) {
                return rootParameter.getParameterValue(0);
            } else {
                return branchRateModel.getBranchRate(tree, tree.getNode(nodeNum));
            }
        }

        @Override
        public void setParameterValue(int dim, double value) {

        }

        @Override
        public void setParameterValueQuietly(int dim, double value) {

        }

        @Override
        public void setParameterValueNotifyChangedAll(int dim, double value) {

        }

        @Override
        public void addBounds(Bounds<Double> bounds) {
        }

        @Override
        public void modelChangedEvent(Model model, Object object, int index) {
            fireParameterChangedEvent();
        }

        @Override
        public void modelRestored(Model model) {

        }
    }

    @Override
    public Mapping getBranchModelMapping(NodeRef branch) {
        return new BranchModel.Mapping() {
            public int[] getOrder() {return new int[] {branch.getNumber()};}

            @Override
            public double[] getWeights() {
                return new double[] {1.0};
            }
        };
    }

    @Override
    public List<SubstitutionModel> getSubstitutionModels() {
        return substitutionModelList;
    }

    public SubstitutionModel getSubstitutionModel(NodeRef branch) {
        return substitutionModelList.get(branch.getNumber());
    }

    @Override
    public SubstitutionModel getRootSubstitutionModel() {
        return substitutionModelList.get(tree.getRoot().getNumber());
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
}
