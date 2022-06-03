/*
 * CountableBranchCategoryProvider.java
 *
 * Copyright (c) 2002-2022 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
import dr.evomodel.tree.TreeModel;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

/**
 * @author Marc A. Suchard
 * @author Philippe Lemey
 */
public class IntegratedTimeDependentModel extends AbstractModel implements ContinuousBranchValueProvider,
        CountableMixtureBranchRates.TimeDependentModel {

    private final TreeModel treeModel;
    private final ParameterPack pack;

    private HashMap<BranchLength, Double> integratedValues;
    private HashMap<BranchLength, Double> savedIntegratedValues;

    public IntegratedTimeDependentModel(TreeModel treeModel, ParameterPack pack) {
        super("integratedBranchValues");
        this.treeModel = treeModel;
        this.pack = pack;

        addModel(treeModel);

        for (Parameter p : pack) {
            addVariable(p);
        }
    }

    private double computeIntegratedValue(double parent, double child) {
        return parent - child;
    }

    @Override
    public double getBranchValue(Tree tree, NodeRef node) {

        double parent = tree.getNodeHeight(tree.getParent(node));
        double child = tree.getNodeHeight(node);
        BranchLength branchLength = new BranchLength(parent, child);

        return integratedValues.computeIfAbsent(branchLength,
                k -> computeIntegratedValue(parent, child));
    }

    @Override
    public double getMidpointValue(Tree tree, NodeRef node, boolean log) {
        double midpoint = getBranchValue(tree, node);

        if (log) {
            return Math.log(midpoint);
        } else {
            return midpoint;
        }
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model != treeModel) {
            throw new IllegalArgumentException("Unknown model");
        }
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if (pack.contains(variable)) {
            integratedValues.clear();
        } else {
            throw new IllegalArgumentException("Unknown variable");
        }
    }

    @Override
    protected void storeState() {
        savedIntegratedValues.clear();
        savedIntegratedValues.putAll(integratedValues);
    }

    @Override
    protected void restoreState() {
        HashMap<BranchLength, Double> tmp = integratedValues;
        integratedValues = savedIntegratedValues;
        savedIntegratedValues = tmp;
    }

    @Override
    protected void acceptState() { }

    static class BranchLength {

        final double parent;
        final double child;

        BranchLength(double parent, double child) {
            this.parent = parent;
            this.child = child;
        }
    }

    public static class ParameterPack implements Iterable<Parameter> {

        final Parameter historicValue;
        final Parameter currentValue;
        final Parameter midTime;
        final Parameter slope;

        final List<Parameter> parameterList = new ArrayList<>();

        public ParameterPack(Parameter historicValue,
                             Parameter currentValue,
                             Parameter midTime,
                             Parameter slope) {
            this.historicValue = historicValue;
            this.currentValue = currentValue;
            this.midTime = midTime;
            this.slope = slope;

            parameterList.add(historicValue);
            parameterList.add(currentValue);
            parameterList.add(midTime);
            parameterList.add(slope);
        }

        public boolean contains(Variable variable) {
            return parameterList.contains((Parameter) variable);
        }

        @Override
        public Iterator<Parameter> iterator() {
            return parameterList.iterator();
        }
    }
}

