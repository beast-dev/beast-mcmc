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
import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.branchratemodel.ArbitraryBranchRates.BranchRateTransform;
import dr.evomodel.branchratemodel.BranchRateModel;

/**
 * @author Marc Suchard
 * @author Xiang Ji
 */
public class BranchParameter extends Parameter.Abstract implements VariableListener, ModelListener {

    final private CompoundParameter parameter;
    final private BranchRateModel branchRateModel;
    final private Parameter rootParameter;
    final private Tree tree;

    public BranchParameter(String name,
                           Tree tree,
                           BranchRateModel branchRateModel,
                           Parameter rootParameter) {
        super(name);

        this.rootParameter = rootParameter;
        this.tree = tree;
        this.branchRateModel = branchRateModel;

        branchRateModel.addModelListener(this);
        this.parameter = constructParameter();

    }

    private CompoundParameter constructParameter() {
        CompoundParameter compoundParameter =  new CompoundParameter(getId() + ".parameter");
        for (int i = 0; i < tree.getNodeCount(); i++) {
            BranchSpecificProxyParameter proxyParameter = new BranchSpecificProxyParameter(branchRateModel, tree, i);
            compoundParameter.addParameter(proxyParameter);
        }
        return compoundParameter;
    }

    public Parameter getRootParameter() {
        return rootParameter;
    }

    public BranchSpecificProxyParameter getParameter(int dim) {
        return (BranchSpecificProxyParameter) parameter.getParameter(dim);
    }

    public BranchRateTransform getTransform() {
        if (branchRateModel instanceof ArbitraryBranchRates) {
            return ((ArbitraryBranchRates) branchRateModel).getTransform();
        } else {
            throw new RuntimeException("Not yet implemented!");
        }
    }

    public Parameter getParameter() {
        return parameter;
    }

    public double[] getParameterValues() {
        return parameter.getParameterValues();
    }

    @Override
    public double getParameterValue(int dim) {

        if (dim == tree.getRoot().getNumber()) {
            return rootParameter.getParameterValue(0);
        } else {
            return branchRateModel.getBranchRate(tree, tree.getNode(dim));
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

    public double getChainGradient(Tree tree, NodeRef node) {
//        final double raw = parameter.getParameterValue(branchRateModel.getParameterIndexFromNode(node));
//        return branchRateModel.getTransform().differential(raw, tree, node);
        throw new RuntimeException("Not yet implemented!");
    }

    @Override
    public void modelChangedEvent(Model model, Object object, int index) {
        fireParameterChangedEvent();
    }

    @Override
    public void modelRestored(Model model) {

    }

    private class BranchSpecificProxyParameter extends Parameter.Proxy {
        private BranchRateModel branchRateModel;
        private final int nodeNum;
        private Tree tree;

        private BranchSpecificProxyParameter(BranchRateModel branchRateModel,
                                             Tree tree,
                                             int nodeNum) {
            super("BranchSpecificProxyParameter." + Integer.toString(nodeNum), 1);
            this.branchRateModel = branchRateModel;
            this.nodeNum = nodeNum;
            this.tree = tree;
        }

        @Override
        public double getParameterValue(int dim) {
            return branchRateModel.getBranchRate(tree, tree.getNode(nodeNum));
        }

        @Override
        public void setParameterValue(int dim, double value) {
//            branchRateModel.setBranchRate(tree, tree.getNode(nodeNum), value);
            throw new RuntimeException("Not yet implemented!");
        }

        @Override
        public void setParameterValueQuietly(int dim, double value) {
            throw new RuntimeException("Not yet implemented!");
        }

        @Override
        public void setParameterValueNotifyChangedAll(int dim, double value) {
            throw new RuntimeException("Not yet implemented!");
        }

        @Override
        public void addBounds(Bounds<Double> bounds) {
//            if (getBounds() == null) {
//                super.addBounds(bounds);
//            }
        }

        @Override
        public Bounds<Double> getBounds() {
//            return branchRateModel.getRateParameter().getBounds();
            return null;
        }
    }
}
