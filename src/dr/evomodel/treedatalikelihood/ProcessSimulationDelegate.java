/*
 * DataSimulationDelegate.java
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
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.continuous.cdi.ContinuousDiffusionIntegrator;
import dr.inference.model.Model;
import dr.inference.model.ModelListener;

import java.util.Arrays;
import java.util.List;

/**
 * ProcessSimulationDelegate - interface for a plugin delegate for data simulation on a tree.
 *
 * @author Andrew Rambaut
 * @author Marc Suchard
 * @version $Id$
 */
public interface ProcessSimulationDelegate extends ProcessOnTreeDelegate {

    void simulate(List<DataLikelihoodDelegate.BranchOperation> branchOperations,
                  List<DataLikelihoodDelegate.NodeOperation> nodeOperations,
                  int rootNodeNumber);

    abstract class AbstractDelegate implements ProcessSimulationDelegate {

        AbstractDelegate(double[] simulation) {
            this.simulation = simulation;
        }

        @Override
        public final TreeTraversal.TraversalType getOptimalTraversalType() {
            return TreeTraversal.TraversalType.PRE_ORDER;
        }

        @Override
        public final void setCallback(TreeDataLikelihood treeDataLikelihood) {
            // TODO Is a call back ever necessary?
        }

        protected final double[] simulation;
    }

    class ConditionalOnTipsDelegate extends AbstractDelegate {

        ConditionalOnTipsDelegate(double[] simulation) {
            super(simulation);
        }

        @Override
        public void simulate(List<BranchOperation> branchOperations, List<NodeOperation> nodeOperations, int rootNodeNumber) {
            Arrays.fill(simulation, 42.0); // TODO Do nothing yet
        }
    }

    class UnconditionalOnTipsDelegate extends AbstractDelegate {

        UnconditionalOnTipsDelegate(double[] simulation) {
            super(simulation);
        }

        @Override
        public void simulate(List<BranchOperation> branchOperations, List<NodeOperation> nodeOperations, int rootNodeNumber) {
            Arrays.fill(simulation, 42.0); // TODO Do nothing yet
        }

    }

    class Test implements ModelListener, TreeTraitProvider { // TODO Rename and move into separate unit

        private final TreeModel treeModel;
        private final String traitName;
        private final TreeTraversal treeTraversalDelegate;
        private final TreeDataLikelihood treeDataLikelihood;
        private final ContinuousDiffusionIntegrator integrator;
        private final ContinuousDataLikelihoodDelegate likelihoodDelegate;
        private final Helper treeTraitProvider = new TreeTraitProvider.Helper();

        private final ConditionalOnTipsDelegate simulationDelegate;

        private final int dimNode;
        private final double[] simulation;

        private boolean validSimulation;

        public Test(String traitName, TreeModel treeModel, TreeDataLikelihood treeDataLikelihood, BranchRateModel branchRateModel,
                    ContinuousDiffusionIntegrator integrator) {
            this(traitName, treeModel, treeDataLikelihood, new TreeTraversal(treeModel, branchRateModel,
                    TreeTraversal.TraversalType.PRE_ORDER), integrator);
        }

        public Test(String traitName, TreeModel treeModel, TreeDataLikelihood treeDataLikelihood, TreeTraversal treeTraversalDelegate,
                    ContinuousDiffusionIntegrator integrator) {

            assert treeDataLikelihood.getDataLikelihoodDelegate() instanceof ContinuousDataLikelihoodDelegate : "Implemented for continuous traits";

            this.traitName = traitName;
            this.treeModel = treeModel;
            this.treeDataLikelihood = treeDataLikelihood;
            this.likelihoodDelegate = (ContinuousDataLikelihoodDelegate) treeDataLikelihood.getDataLikelihoodDelegate();
            this.integrator = integrator;
            this.treeTraversalDelegate = treeTraversalDelegate;

            treeDataLikelihood.addModelListener(this);

            dimNode = likelihoodDelegate.getNumTraits() * likelihoodDelegate.getDimTrait();
            simulation = new double[treeModel.getNodeCount() * dimNode];

            simulationDelegate = new ConditionalOnTipsDelegate(simulation);

            addTraitToProvider();

            validSimulation = false;
        }

        private void addTraitToProvider() {
            if (treeTraitProvider.getTreeTrait(traitName) == null) {
                treeTraitProvider.addTrait(
                        new TreeTrait.DA() {
                            public String getTraitName() {
                                return traitName;
                            }

                            public Intent getIntent() {
                                return Intent.NODE;
                            }

                            public Class getTraitClass() {
                                return Double.class;
                            }

                            public double[] getTrait(Tree tree, NodeRef node) {
                                assert tree == treeModel : "Bad tree";

                                return getTraitForNode(tree, node);
                            }
                        });
            }
        }

        private double[] getTraitForAllTips() {

            cacheSimulatedTraits(null);

            final int length = dimNode * treeTraversalDelegate.getTree().getExternalNodeCount();
            double[] trait = new double[length];
            System.arraycopy(simulation, 0, trait, 0, length);

            return trait;
        }

        private double[] getTraitForNode(final Tree tree, final NodeRef node) {

            cacheSimulatedTraits(null);

            double[] trait = new double[dimNode];
            System.arraycopy(simulation, node.getNumber() * dimNode, trait, 0, dimNode);

            return trait;
        }

        private void cacheSimulatedTraits(final NodeRef node) {
            treeDataLikelihood.getLogLikelihood(); // Ensure likelihood is up-to-date

            if (!validSimulation) {
                simulateTraits(node);
                validSimulation = true;
            }
        }

        private void simulateTraits(final NodeRef targetNode) {

            if (targetNode == null) {
                treeTraversalDelegate.updateAllNodes(); // TODO depends on targetNode
            } else {
                throw new RuntimeException("Not yet implemented");
            }

            treeTraversalDelegate.dispatchTreeTraversalCollectBranchAndNodeOperations();

            List<BranchOperation> branchOperations = treeTraversalDelegate.getBranchOperations();
            List<NodeOperation> nodeOperations = treeTraversalDelegate.getNodeOperations();

            final NodeRef root = treeModel.getRoot();
            simulationDelegate.simulate(branchOperations, nodeOperations, root.getNumber());

            treeTraversalDelegate.setAllNodesUpdated();
        }

        @Override
        public TreeTrait[] getTreeTraits() {
            return treeTraitProvider.getTreeTraits();
        }

        @Override
        public TreeTrait getTreeTrait(String key) {
            return treeTraitProvider.getTreeTrait(key);
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
}
