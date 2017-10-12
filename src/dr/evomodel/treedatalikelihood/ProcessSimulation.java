/*
 * ProcessSimulation.java
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

package dr.evomodel.treedatalikelihood;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evolution.tree.TreeTraitProvider;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.treedatalikelihood.continuous.cdi.ContinuousDiffusionIntegrator;
import dr.evomodel.treedatalikelihood.preorder.ProcessSimulationDelegate;
import dr.inference.model.Model;
import dr.inference.model.ModelListener;

import java.util.List;

/**
 * @author Marc A. Suchard
 * @author Andrew Rambaut
 * @author Philippe Lemey
 */
public class ProcessSimulation implements ModelListener, TreeTraitProvider {

    private final Tree tree;
    private final String name;

    private final SimulationTreeTraversal treeTraversalDelegate;
    private final TreeDataLikelihood treeDataLikelihood;
    private final DataLikelihoodDelegate likelihoodDelegate;
    private final ProcessSimulationDelegate simulationDelegate;

    private final int[] operations;

    private boolean validSimulation;

    public ProcessSimulation(String name, TreeDataLikelihood treeDataLikelihood,
                             ProcessSimulationDelegate simulationDelegate) {

        this.name = name;
        this.treeDataLikelihood = treeDataLikelihood;
        this.tree = treeDataLikelihood.getTree();

        BranchRateModel branchRateModel = treeDataLikelihood.getBranchRateModel();
        treeTraversalDelegate = new SimulationTreeTraversal(tree, branchRateModel,
                simulationDelegate.getOptimalTraversalType());

        this.likelihoodDelegate = treeDataLikelihood.getDataLikelihoodDelegate();
        treeDataLikelihood.addModelListener(this);

        this.simulationDelegate = simulationDelegate;
        simulationDelegate.setCallback(this);

        this.operations = new int[tree.getNodeCount() * ContinuousDiffusionIntegrator.OPERATION_TUPLE_SIZE];

        validSimulation = false;
    }

    public final void cacheSimulatedTraits(final NodeRef node) {

        treeDataLikelihood.getLogLikelihood(); // Ensure likelihood is up-to-date

        if (!validSimulation) {
            simulateTraits(node);
            validSimulation = true;
        }
    }

    private final void simulateTraits(final NodeRef targetNode) {

        if (targetNode == null) {
            treeTraversalDelegate.updateAllNodes();
        } else {
            treeTraversalDelegate.updateAllNodes(); // TODO Fix - depends on targetNode
        }

        treeTraversalDelegate.dispatchTreeTraversalCollectBranchAndNodeOperations();
        final NodeRef root = tree.getRoot();

        boolean NEW = false;

        if (NEW) {
            int count = convertOperationsToList(treeTraversalDelegate.getNodeOperations(), operations);
            simulationDelegate.simulate(operations, count, root.getNumber());
        } else {
            simulationDelegate.simulate(treeTraversalDelegate, root.getNumber());
        }

        treeTraversalDelegate.setAllNodesUpdated();
    }

    private int convertOperationsToList(final List<ProcessOnTreeDelegate.NodeOperation> nodeOperations,
                                         int[] operations) {

        int k = 0;
//        for (ProcessOnTreeDelegate.NodeOperation op : nodeOperations) {
//            operations[k    ] = likelihoodDelegate.getActiveNodeIndex(op.getNodeNumber());
//            operations[k + 1] = likelihoodDelegate.getActiveNodeIndex(op.getLeftChild());    // source node 1
//            operations[k + 2] = likelihoodDelegate.getActiveMatrixIndex(op.getLeftChild());  // source matrix 1
//            operations[k + 3] = likelihoodDelegate.getActiveNodeIndex(op.getRightChild());   // source node 2
//            operations[k + 4] = likelihoodDelegate.getActiveMatrixIndex(op.getRightChild()); // source matrix 2
//
//            k += ContinuousDiffusionIntegrator.OPERATION_TUPLE_SIZE;
//        }

        return nodeOperations.size();
    }

    @Override
    public TreeTrait[] getTreeTraits() {
        return simulationDelegate.getTreeTraits();
    }

    @Override
    public TreeTrait getTreeTrait(String key) {
        return simulationDelegate.getTreeTrait(key);
    }

    @Override
    public void modelChangedEvent(Model model, Object object, int index) {
        assert model == treeDataLikelihood : "Invalid model";

        validSimulation = false;
    }

    @Override
    public void modelRestored(Model model) {
        // Do nothing
    }
}
