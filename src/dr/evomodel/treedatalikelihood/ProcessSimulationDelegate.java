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

import dr.evolution.tree.*;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.treedatalikelihood.continuous.ConjugateRootTraitPrior;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.continuous.ContinuousRateTransformation;
import dr.evomodel.treedatalikelihood.continuous.ContinuousTraitDataModel;

import java.util.Arrays;
import java.util.List;

/**
 * ProcessSimulationDelegate - interface for a plugin delegate for data simulation on a tree.
 *
 * @author Andrew Rambaut
 * @author Marc Suchard
 * @version $Id$
 */
public interface ProcessSimulationDelegate extends ProcessOnTreeDelegate, TreeTraitProvider {

    void simulate(List<DataLikelihoodDelegate.BranchOperation> branchOperations,
                  List<DataLikelihoodDelegate.NodeOperation> nodeOperations,
                  int rootNodeNumber);

    void setCallback(ProcessSimulation simulationProcess);

    abstract class AbstractDelegate implements ProcessSimulationDelegate {

        AbstractDelegate(String name, Tree tree) {
            this.name = name;
            this.tree = tree;
            constructTraits(treeTraitHelper);
        }

        abstract void constructTraits(Helper treeTraitHelper);

        @Override
        public final TreeTraversal.TraversalType getOptimalTraversalType() {
            return TreeTraversal.TraversalType.PRE_ORDER;
        }

        @Override
        public final void setCallback(ProcessSimulation simulationProcess) {
            this.simulationProcess = simulationProcess;
        }

        @Override
        public final TreeTrait[] getTreeTraits() {
            return treeTraitHelper.getTreeTraits();
        }

        @Override
        public final TreeTrait getTreeTrait(String key) {
            return treeTraitHelper.getTreeTrait(key);
        }

        protected final TreeTraitProvider.Helper treeTraitHelper = new Helper();

        protected ProcessSimulation simulationProcess = null;

        protected final Tree tree;
        protected final String name;

    }

    abstract class AbstractContinuousTraitDelegate extends AbstractDelegate {

        protected final int dimTrait;
        protected final int numTraits;
        protected final int dimNode;


        AbstractContinuousTraitDelegate(String name,
                                        MultivariateTraitTree tree,
                                        MultivariateDiffusionModel diffusionModel,
                                        ContinuousTraitDataModel dataModel,
                                        ConjugateRootTraitPrior rootPrior,
                                        ContinuousRateTransformation rateTransformation,
                                        BranchRateModel rateModel,
                                        ContinuousDataLikelihoodDelegate likelihoodDelegate) {
            super(name, tree);

            dimTrait = likelihoodDelegate.getTraitDim();
            numTraits = likelihoodDelegate.getTraitCount();
            dimNode = dimTrait * numTraits;

            this.buffer = new double[dimNode * tree.getNodeCount()];
        }

        @Override
        void constructTraits(Helper treeTraitHelper) {

            TreeTrait.DA baseTrait = new TreeTrait.DA() {

                public String getTraitName() {
                    return name;
                }

                public Intent getIntent() {
                    return Intent.NODE;
                }

                public double[] getTrait(Tree t, NodeRef node) {

                    assert t == tree;
                    return getTraitForNode(node);
                }
            };

            treeTraitHelper.addTrait(baseTrait);
        }

        private double[] getTraitForAllTips() { // TODO To be used as a GaussianProcessRandomGenerator

            assert simulationProcess != null;

            simulationProcess.cacheSimulatedTraits(null);

            final int length = dimNode * tree.getExternalNodeCount();
            double[] trait = new double[length];
            System.arraycopy(buffer, 0, trait, 0, length);

            return trait;
        }

        private double[] getTraitForNode(final NodeRef node) {

            assert simulationProcess != null;

            simulationProcess.cacheSimulatedTraits(null);

            double[] trait = new double[dimNode];
            System.arraycopy(buffer, node.getNumber() * dimNode, trait, 0, dimNode);

            return trait;
        }

        protected final double[] buffer;
    }

    class ConditionalOnTipsDelegate extends AbstractContinuousTraitDelegate {

        public ConditionalOnTipsDelegate(String name,
                                         MultivariateTraitTree tree,
                                         MultivariateDiffusionModel diffusionModel,
                                         ContinuousTraitDataModel dataModel,
                                         ConjugateRootTraitPrior rootPrior,
                                         ContinuousRateTransformation rateTransformation,
                                         BranchRateModel rateModel,
                                         ContinuousDataLikelihoodDelegate likelihoodDelegate) {
            super(name, tree, diffusionModel, dataModel, rootPrior, rateTransformation, rateModel, likelihoodDelegate);
        }

        @Override
        public void simulate(List<BranchOperation> branchOperations, List<NodeOperation> nodeOperations, int rootNodeNumber) {
            Arrays.fill(buffer, 42.0); // TODO Do nothing yet
        }
    }

    class UnconditionalOnTipsDelegate extends AbstractContinuousTraitDelegate {

        public UnconditionalOnTipsDelegate(String name,
                                           MultivariateTraitTree tree,
                                           MultivariateDiffusionModel diffusionModel,
                                           ContinuousTraitDataModel dataModel,
                                           ConjugateRootTraitPrior rootPrior,
                                           ContinuousRateTransformation rateTransformation,
                                           BranchRateModel rateModel,
                                           ContinuousDataLikelihoodDelegate likelihoodDelegate) {
            super(name, tree, diffusionModel, dataModel, rootPrior, rateTransformation, rateModel, likelihoodDelegate);
        }

        @Override
        public void simulate(List<BranchOperation> branchOperations, List<NodeOperation> nodeOperations, int rootNodeNumber) {
            Arrays.fill(buffer, 42.0); // TODO Do nothing yet
        }

    }

}
