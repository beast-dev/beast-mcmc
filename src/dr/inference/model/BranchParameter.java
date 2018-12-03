/*
 * BranchParameter.java
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

package dr.inference.model;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.branchratemodel.ArbitraryBranchRates.BranchRateTransform;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.TreeParameterModel;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Marc Suchard
 * @author Xiang Ji
 */
public class BranchParameter extends Parameter.Abstract implements VariableListener, ModelListener {

    final private CompoundParameter parameter;
    final private Parameter rootParameter;
    final private BranchRateTransform transform;
    final private TreeModel tree;
    final private TreeParameterModel indexHelper;
    private Boolean addedTransformedParameter;
    private List<IndividualBranchParameter> transformedParameter = new ArrayList<IndividualBranchParameter>();

    public BranchParameter(String name,
                           CompoundParameter parameter,
                           Parameter rootParameter,
                           TreeModel tree,
                           BranchRateTransform transform) {

        super(name);

        this.parameter =  parameter;
        this.rootParameter = rootParameter;
        this.transform = transform;
        this.tree = tree;
        this.indexHelper = new TreeParameterModel(tree, parameter, false, TreeTrait.Intent.BRANCH);
        this.parameter.addVariableListener(this);
        this.rootParameter.addVariableListener(this);
        if (transform instanceof Model) {
            ((Model) transform).addModelListener(this);
        }

        this.addedTransformedParameter = false;
    }

    public void addTransformedParameterList(List<IndividualBranchParameter> transformedParameter) {
        if (addedTransformedParameter) {
            throw new RuntimeException("Should be called only once.");
        } else {
            if (transformedParameter.size() != tree.getNodeCount()) {
                throw new RuntimeException("Size mismatch!");
            }
            this.transformedParameter = transformedParameter;
            addedTransformedParameter = true;
        }
    }

    public IndividualBranchParameter getParameter(int dim) {
        return transformedParameter.get(dim);
    }

    public BranchRateTransform getTransform() {
        return transform;
    }

    public Parameter getParameter() {
        return parameter;
    }

    public double[] getParameterValues() {
        return getBranchParameterValues();
    }

    @Override
    public double getParameterValue(int dim) {

        if (dim == tree.getNodeCount() - 1) {
            return rootParameter.getParameterValue(0);
        } else {
            return transform.transform(parameter.getParameterValue(dim), tree, tree.getNode(indexHelper.getNodeNumberFromParameterIndex(dim)));
        }

    }

    @Override
    public void setParameterValue(int dim, double value) {
        parameter.setParameterValue(dim, value);
    }

    @Override
    public void setParameterValueQuietly(int dim, double value) {
        parameter.setParameterValueQuietly(dim, value);
    }

    @Override
    public void setParameterValueNotifyChangedAll(int dim, double value) {
        parameter.setParameterValueNotifyChangedAll(dim, value);
    }

    @Override
    public String getParameterName() {
        String name = getId();
        if (name == null) {
            name = "BranchParameter." + parameter.getParameterName();
        }
        return name;
    }

    @Override
    public void addBounds(Bounds<Double> bounds) {
        parameter.addBounds(bounds);
    }

    @Override
    public Bounds<Double> getBounds() {
        return parameter.getBounds();
    }

    @Override
    public void addDimension(int index, double value) {
        throw new RuntimeException("Dimension should not be changed.");
    }

    @Override
    public double removeDimension(int index) {
        throw new RuntimeException("Dimension should not be changed.");
    }

    @Override
    protected void storeValues() {
        parameter.storeParameterValues();
    }

    @Override
    protected void restoreValues() {
        parameter.restoreParameterValues();
    }

    @Override
    protected void acceptValues() {
        parameter.acceptParameterValues();
    }

    @Override
    protected void adoptValues(Parameter source) {
        parameter.adoptParameterValues(source);
    }

    @Override
    public void variableChangedEvent(Variable variable, int index, ChangeType type) {
        fireParameterChangedEvent(index, type);
    }

    public int getDimension() {
        return tree.getNodeCount() - 1;
    }

    public double[] getBranchParameterValues() {
        double[] copyOfValues = new double[tree.getNodeCount() - 1];
        for (int i = 0; i < copyOfValues.length; i++) {
            copyOfValues[i] = getParameterValue(indexHelper.getNodeNumberFromParameterIndex(i));
        }
        return copyOfValues;
    }

    public double getChainGradient(Tree tree, NodeRef node) {
        final double raw = parameter.getParameterValue(node.getNumber());
        return transform.differential(raw, tree, node);
    }

    @Override
    public void modelChangedEvent(Model model, Object object, int index) {
        fireParameterChangedEvent();
    }

    @Override
    public void modelRestored(Model model) {

    }

    public static class IndividualBranchParameter extends Parameter.Abstract implements VariableListener, ModelListener {

        final private Parameter parameter;
        final private BranchParameter branchParameter;
        final private int nodeNum;

        public IndividualBranchParameter(BranchParameter branchParameter, int nodeNum, Parameter parameter) {
            this.branchParameter = branchParameter;
            this.nodeNum = nodeNum;
            this.parameter = parameter;
            if (!(parameter.getDimension() == 1)) {
                throw new RuntimeException("Individual parameter can only be one dimensional.");
            }
            this.parameter.addVariableListener(this);
            this.branchParameter.addParameterListener(this);

        }

        public BranchParameter getBranchParameter() {
            return branchParameter;
        }

        @Override
        public double getParameterValue(int dim) {
            if (dim > 0) {
                throw new RuntimeException("Should be one dimensional!");
            }
            return branchParameter.getParameterValue(nodeNum);
        }

        @Override
        public void setParameterValue(int dim, double value) {
            branchParameter.setParameterValue(nodeNum, value);
        }

        @Override
        public void setParameterValueQuietly(int dim, double value) {
            branchParameter.setParameterValueQuietly(nodeNum, value);
        }

        @Override
        public void setParameterValueNotifyChangedAll(int dim, double value) {
            branchParameter.setParameterValueNotifyChangedAll(nodeNum, value);
        }

        @Override
        public String getParameterName() {
            return branchParameter.getParameterName() + "." + nodeNum;
        }

        @Override
        public void addBounds(Bounds<Double> bounds) {
            parameter.addBounds(bounds);
        }

        @Override
        public Bounds<Double> getBounds() {
            return parameter.getBounds();
        }

        @Override
        public void addDimension(int index, double value) {
            throw new RuntimeException("Fixed dimension should not be changed.");
        }

        @Override
        public double removeDimension(int index) {
            throw new RuntimeException("Fixed dimension should not be changed.");
        }

        @Override
        protected void storeValues() {
            parameter.storeParameterValues();
        }

        @Override
        protected void restoreValues() {
            parameter.restoreParameterValues();
        }

        @Override
        protected void acceptValues() {
            parameter.acceptParameterValues();
        }

        @Override
        protected void adoptValues(Parameter source) {
            parameter.adoptParameterValues(source);
        }

        @Override
        public void variableChangedEvent(Variable variable, int index, ChangeType type) {
            fireParameterChangedEvent(index, type);
        }

        @Override
        public void modelChangedEvent(Model model, Object object, int index) {
            fireParameterChangedEvent();
        }

        @Override
        public void modelRestored(Model model) {

        }
    }
}
