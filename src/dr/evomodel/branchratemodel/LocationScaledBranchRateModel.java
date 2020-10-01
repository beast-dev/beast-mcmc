/*
 * RandomLocalClockModel.java
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
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.util.Citable;
import dr.util.Citation;

import java.util.ArrayList;
import java.util.List;
import java.util.function.DoubleBinaryOperator;

/**
 * @author Marc A. Suchard
 * @author Alexander Fisher
 */

public class LocationScaledBranchRateModel extends AbstractBranchRateModel
        implements DifferentiableBranchRates, Citable {

    private final TreeModel treeModel;
    private final BranchRateModel branchRateModel;
    private final DifferentiableBranchRates differentiableBranchRateModel;
    private final BranchSpecificFixedEffects location;

    public LocationScaledBranchRateModel(TreeModel treeModel,
                                         BranchRateModel branchRateModel,
                                         BranchSpecificFixedEffects location) {

        super("name");

        this.treeModel = treeModel;
        this.branchRateModel = branchRateModel;
        this.differentiableBranchRateModel = (branchRateModel instanceof DifferentiableBranchRates) ?
                (DifferentiableBranchRates) branchRateModel : null;
        this.location = location;

        addModel(treeModel);
        addModel(branchRateModel);
        if (location instanceof Model) {
            addModel((Model) location);
        }
    }

    public void handleModelChangedEvent(Model model, Object object, int index) {
        fireModelChanged();
    }

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        fireModelChanged();
    }

    protected void storeState() { }

    protected void restoreState() { }

    protected void acceptState() { }

    @Override
    public double getBranchRateDifferential(final Tree tree, final NodeRef node) {
        checkDifferentiability();
        return location.getEffect(tree, node) * differentiableBranchRateModel.getBranchRateDifferential(tree, node);
    }

    @Override
    public double getBranchRateSecondDifferential(Tree tree, NodeRef node) {
        checkDifferentiability();
        return location.getEffect(tree, node) * differentiableBranchRateModel.getBranchRateSecondDifferential(tree, node);
    }

    @Override
    public Parameter getRateParameter() {
        checkDifferentiability();
        return differentiableBranchRateModel.getRateParameter();
    }

    public Tree getTree() {
        return treeModel;
    }

    @Override
    public int getParameterIndexFromNode(NodeRef node) {
        checkDifferentiability();
        return differentiableBranchRateModel.getParameterIndexFromNode(node);
    }

    private void checkDifferentiability() {
        if (differentiableBranchRateModel == null) {
            throw new RuntimeException("Non-differentiable base BranchRateModel");
        }
    }

    @Override
    public ArbitraryBranchRates.BranchRateTransform getTransform() {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double[] updateGradientLogDensity(double[] gradient, double[] value, int from, int to) {
//        throw new RuntimeException("Not yet implemented");
        return gradient;
    }

    @Override
    public double[] updateDiagonalHessianLogDensity(double[] diagonalHessian, double[] gradient, double[] value, int from, int to) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double mapReduceOverRates(NodeRateMap map, DoubleBinaryOperator reduce, double initial) {
        checkDifferentiability();
        return differentiableBranchRateModel.mapReduceOverRates(map, reduce, initial);
    }

    @Override
    public void forEachOverRates(NodeRateMap map) {
        checkDifferentiability();
        differentiableBranchRateModel.forEachOverRates(map);
    }

    @Override
    public double getBranchRate(final Tree tree, final NodeRef node) {
        assert tree == treeModel;
        return location.getEffect(tree, node) * branchRateModel.getBranchRate(tree, node);
    }

    public double getUntransformedBranchRate(Tree tree, NodeRef node) {
        // returns the rate scaled by the location
//        return getBranchRate(tree, node);

        //returns just the rate
        return branchRateModel.getBranchRate(tree, node);
    }

    public double getPriorRateAsIncrement(Tree tree){
        return Math.log(location.getEffect(tree, null));
    }

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.MOLECULAR_CLOCK;
    }

    @Override
    public String getDescription() {
        String description =
                (branchRateModel instanceof Citable) ?
                        ((Citable) branchRateModel).getDescription() :
                        "Unknown clock model";

        description += " with scaling-by-tree-time";
        return description;
    }

    @Override
    public List<Citation> getCitations() {
        List<Citation> list = 
                (branchRateModel instanceof Citable) ?
                        new ArrayList<>(((Citable) branchRateModel).getCitations()) :
                        new ArrayList<>();
        // TODO
        return list;
    }
}