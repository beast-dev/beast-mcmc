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
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.TreeParameterModel;
import dr.evomodel.treedatalikelihood.SimulationTreeTraversal;
import dr.evomodel.treedatalikelihood.preorder.LocalBranchRateDelegate;
import dr.evomodelxml.branchratemodel.LocalBranchRatesParser;
import dr.inference.model.Model;
import dr.inference.model.Parameter;

/**
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public class LocalBranchRates extends ArbitraryBranchRates {

    private TreeModel tree;
    private TreeParameterModel ratesMultiplier;
    private Parameter multiplierParameter;

    private TreeParameterModel branchRates;
    private double normalizingConstant;
    private LocalBranchRateDelegate branchRateDelegate;

    private final BranchRateTransform transform;
    private final SimulationTreeTraversal treeTraversalDelegate;
    private final int[] operations;


    public LocalBranchRates(TreeModel tree,
                            Parameter multiplierParameter,
                            BranchRateTransform transform) {

        super(tree, multiplierParameter, new BranchRateTransform.None(), false);
        this.tree = tree;
        this.multiplierParameter = multiplierParameter;
        this.ratesMultiplier = new TreeParameterModel(tree, this.multiplierParameter, false, Intent.BRANCH);
        this.branchRates = new TreeParameterModel(tree, new Parameter.Default(tree.getNodeCount()), true, Intent.NODE);
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
        double totalBranchLength = 0.0;
        double totalRateBranchLengthProduct = 0.0;
        for (int i = 0; i < tree.getNodeCount(); i++) {
            NodeRef node = tree.getNode(i);
            if (!tree.isRoot(node)) {
                totalBranchLength += tree.getBranchLength(node);
                totalRateBranchLengthProduct += tree.getBranchLength(node) * ratesMultiplier.getNodeValue(tree, node);
            }
        }
        normalizingConstant = totalBranchLength / totalRateBranchLengthProduct;
    }

    private void updateBranchRates() {
        simulateTraits(null);
        updateNormalizingConstant();
    }

    public void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == ratesMultiplier) {
            fireModelChanged(object, index);
        }
    }
}
