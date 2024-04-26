/*
 * LocalBranchRates.java
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
import dr.evomodel.tree.TreeParameterModel;
import dr.evomodel.treedatalikelihood.SimulationTreeTraversal;
import dr.evomodel.treedatalikelihood.preorder.LocalBranchRateDelegate;
import dr.evomodelxml.branchratemodel.LocalBranchRatesParser;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.xml.Reportable;

/**
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
@Deprecated   // What is this class?
public class LocalBranchRates extends ArbitraryBranchRates implements Reportable {

    private TreeModel tree;
    private TreeParameterModel ratesMultiplier;
    private Parameter multiplierParameter;

    private TreeParameterModel branchRates;
    private Parameter branchRateParameter;
    private double normalizingConstant;
    private LocalBranchRateDelegate branchRateDelegate;

    private final BranchRateTransform transform;
    private final SimulationTreeTraversal treeTraversalDelegate;
    private final int[] operations;


    public LocalBranchRates(TreeModel tree,
                            Parameter multiplierParameter,
                            BranchRateTransform transform) {
        this(LocalBranchRatesParser.SHRINKAGE_BRANCH_RATES, tree, multiplierParameter, transform);
    }

    public LocalBranchRates(String name,
                            TreeModel tree,
                            Parameter multiplierParameter,
                            BranchRateTransform transform) {

        super(name, tree, multiplierParameter, new BranchRateTransform.None(), false);
        this.tree = tree;
        this.multiplierParameter = multiplierParameter;
        this.ratesMultiplier = new TreeParameterModel(tree, this.multiplierParameter, false, Intent.BRANCH);
        this.branchRateParameter = new Parameter.Default(tree.getNodeCount());
        this.branchRates = new TreeParameterModel(tree, branchRateParameter, true, Intent.NODE);
        this.branchRateDelegate = new LocalBranchRateDelegate(LocalBranchRatesParser.SHRINKAGE_BRANCH_RATES,
                ratesMultiplier, branchRates);
        this.transform = transform;
        treeTraversalDelegate = new SimulationTreeTraversal(tree, this, branchRateDelegate.getOptimalTraversalType());
        this.operations = new int[tree.getNodeCount() * branchRateDelegate.getSingleOperationSize()];
        updateBranchRates();

        addModel(ratesMultiplier);

        if (transform instanceof Model) {
            addModel((Model) transform);
        }

        if (!(transform instanceof BranchRateTransform.None)) {
            throw new RuntimeException("Only tested without transformation.");
        }
    }

    private void simulateTraits(final NodeRef targetNode) {
        if (targetNode == null) {
            treeTraversalDelegate.updateAllNodes();
        } else {
            treeTraversalDelegate.updateAllNodes(); // TODO: Fix - depends on targetNode as in ProcessSimulation
        }

        treeTraversalDelegate.dispatchTreeTraversalCollectBranchAndNodeOperations();
        int count = branchRateDelegate.vectorizeNodeOperations(treeTraversalDelegate.getNodeOperations(),operations);

        final NodeRef root = tree.getRoot();
        branchRateDelegate.simulate(operations, count, root.getNumber());

        treeTraversalDelegate.setAllNodesUpdated();
    }

    private void updateNormalizingConstant() {
        double totalBranchLength = getTotalBranchLength();
        double totalRateBranchLengthProduct = getTotalRateBranchLengthProduct();

        normalizingConstant = totalBranchLength / totalRateBranchLengthProduct;
//        for (int i = 0; i < tree.getNodeCount(); i++) {
//            final double normalizedRate = branchRateParameter.getParameterValue(i) * normalizingConstant;
//            branchRateParameter.setParameterValue(i, normalizedRate);
//        }
    }

    private double getTotalBranchLength() {
        double totalBranchLength = 0.0;
        for (int i = 0; i < tree.getNodeCount(); i++) {
            NodeRef node = tree.getNode(i);
            if (!tree.isRoot(node)) {
                totalBranchLength += tree.getBranchLength(node);
            }
        }
        return totalBranchLength;
    }

    private double getTotalRateBranchLengthProduct() {
        double totalRateBranchLengthProduct = 0.0;
        for (int i = 0; i < tree.getNodeCount(); i++) {
            NodeRef node = tree.getNode(i);
            if (!tree.isRoot(node)) {
                totalRateBranchLengthProduct += tree.getBranchLength(node) * ratesMultiplier.getNodeValue(tree, node);
            }
        }
        return totalRateBranchLengthProduct;
    }

    private void updateBranchRates() {
        simulateTraits(null);
        updateNormalizingConstant();
    }

    public void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == ratesMultiplier) {
            fireModelChanged();
            updateBranchRates();
        }
    }

    @Override
    public double getBranchRate(final Tree tree, final NodeRef node) {
        return transform.transform(branchRates.getNodeValue(tree, node), tree, node);
    }

    @Override
    public void setBranchRate(Tree tree, NodeRef node, double value) {
        throw new RuntimeException("This function should not be called.");
    }

    @Override
    public String getReport() {
        StringBuilder sb = new StringBuilder();
        sb.append("Local branch rates:\n").append(new dr.math.matrixAlgebra.Vector(branchRateParameter.getParameterValues()));
        sb.append("\n");
        return sb.toString();
    }

    @Override
    public double getBranchRateDifferential(Tree tree, NodeRef node) {
        final double multiplier = ratesMultiplier.getNodeValue(tree, node);
        return branchRates.getNodeValue(tree, node) / multiplier;
    }

    public double getBranchRateSecondDifferential(Tree tree, NodeRef node) {
        throw new RuntimeException("Not yet implemented.");
    }
}
