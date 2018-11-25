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

import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.ArbitraryBranchRates.BranchRateTransform;
import dr.evomodel.tree.TreeParameterModel;


public class BranchParameter extends Parameter.Abstract implements VariableListener {

    final CompoundParameter parameter;
    final BranchRateTransform transform;
    final Tree tree;
    final TreeParameterModel indexHelper;

    public BranchParameter(CompoundParameter parameter, Tree tree, BranchRateTransform transform, TreeParameterModel indexHelper) {

        this.parameter = parameter;
        this.transform = transform;
        this.tree = tree;
        this.indexHelper = indexHelper;
        this.parameter.addVariableListener(this);

    }

    @Override
    public double getParameterValue(int dim) {
        return transform.transform(parameter.getParameterValue(dim), tree,
                tree.getNode(indexHelper.getNodeNumberFromParameterIndex(dim)));
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

    public static class IndividualBranchParameter extends Parameter.Abstract implements VariableListener {

        final private Parameter parameter;
        final private BranchParameter branchParameter;
        final private int parameterIndex;

        public IndividualBranchParameter(BranchParameter branchParameter, int parameterIndex, Parameter parameter) {
            this.branchParameter = branchParameter;
            this.parameterIndex = parameterIndex;
            this.parameter = parameter;
            if (!(parameter.getDimension() == 1)) {
                throw new RuntimeException("Individual parameter can only be one dimensional.");
            }
            this.parameter.addVariableListener(this);
        }

        @Override
        public double getParameterValue(int dim) {
            if (dim > 0) {
                throw new RuntimeException("Should be one dimensional!");
            }
            return branchParameter.getParameterValue(parameterIndex);
        }

        @Override
        public void setParameterValue(int dim, double value) {
            branchParameter.setParameterValue(parameterIndex, value);
        }

        @Override
        public void setParameterValueQuietly(int dim, double value) {
            branchParameter.setParameterValueQuietly(parameterIndex, value);
        }

        @Override
        public void setParameterValueNotifyChangedAll(int dim, double value) {
            branchParameter.setParameterValueNotifyChangedAll(parameterIndex, value);
        }

        @Override
        public String getParameterName() {
            return branchParameter.getParameterName() + "." + parameterIndex;
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
    }
}
