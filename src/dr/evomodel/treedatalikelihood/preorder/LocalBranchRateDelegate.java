/*
 * LocalBranchRateDelegate.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.treedatalikelihood.preorder;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.tree.TreeParameterModel;
import dr.inference.model.Model;

import java.util.List;

import static dr.evolution.tree.TreeTrait.DA.factory;


/**
 * @author Xiang Ji
 * @author Marc Suchard
 */
@Deprecated // What is this class?
public class LocalBranchRateDelegate extends ProcessSimulationDelegate.AbstractDelegate
        implements TreeTrait.TraitInfo<double[]> {

    public static String PRODUCT_BRANCH_RATE_TRAIT_NAME="ProductBranchRate";

    private final TreeParameterModel rateMultipliers;
    private final TreeParameterModel branchRates;


    public LocalBranchRateDelegate(String name,
                                   TreeParameterModel rateMultipliers,
                                   TreeParameterModel branchRates) {
        super(name, rateMultipliers.getTreeModel());
        this.rateMultipliers = rateMultipliers;
        this.branchRates = branchRates;
    }

    @Override
    public String getTraitName() {
        return PRODUCT_BRANCH_RATE_TRAIT_NAME;
    }

    @Override
    public TreeTrait.Intent getTraitIntent() {
        return TreeTrait.Intent.BRANCH;
    }

    @Override
    public Class getTraitClass() {
        return double[].class;
    }

    @Override
    public double[] getTrait(Tree tree, NodeRef node) {
        final int branchNumber = tree.getNodeCount() - 1;
        double[] trait = new double[branchNumber];
        for (int i = 0; i < branchNumber; i++) {
            trait[i] = branchRates.getNodeValue(tree, tree.getNode(rateMultipliers.getNodeNumberFromParameterIndex(i)));
        }
        return trait;
    }

    @Override
    public boolean isTraitLoggable() {
        return false;
    }

    @Override
    public int vectorizeNodeOperations(List<NodeOperation> nodeOperations, int[] operations) {
        int k = 0;
        for (NodeOperation tmpNodeOperation : nodeOperations) {
            operations[k++] = tmpNodeOperation.getLeftChild();
        }
        return nodeOperations.size();
    }

    @Override
    public int getSingleOperationSize() {
        return 1;
    }

    @Override
    protected void constructTraits(Helper treeTraitHelper) {
        treeTraitHelper.addTrait(factory(this));
    }

    @Override
    public void simulate(final int[] operations, final int operationCount, final int rootNodeNumber) {
        setupStatistics();

        simulateRoot(rootNodeNumber);

        int k = 0;
        for (int i = 0; i < operationCount; ++i) {
            simulateNode(
                    operations[k    ]
            );
            k += getSingleOperationSize();

        }
    }

    private void simulateNode(int nodeNum) {
        final NodeRef currentNode = tree.getNode(nodeNum);
        final NodeRef parentNode = tree.getParent(currentNode);
        final double parentRate = branchRates.getNodeValue(tree, parentNode);
        final double currentRate = parentRate * rateMultipliers.getNodeValue(tree, tree.getNode(nodeNum));
        branchRates.setNodeValue(tree, currentNode, currentRate);
    }

    @Override
    protected void setupStatistics() {

    }

    @Override
    protected void simulateRoot(int rootNumber) {
        branchRates.setNodeValue(tree, tree.getNode(rootNumber), 1.0);
    }

    @Override
    protected void simulateNode(int v0, int v1, int v2, int v3, int v4) {
        throw new RuntimeException("This function should not be called!");
    }

    @Override
    public void modelChangedEvent(Model model, Object object, int index) {

    }

    @Override
    public void modelRestored(Model model) {

    }
}
