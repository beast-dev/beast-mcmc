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
import dr.evomodel.treedatalikelihood.continuous.cdi.ContinuousDiffusionIntegrator;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.math.matrixAlgebra.CholeskyDecomposition;
import dr.math.matrixAlgebra.IllegalDimension;
import dr.math.matrixAlgebra.Matrix;
import dr.math.matrixAlgebra.SymmetricMatrix;

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

    void simulate(List<DataLikelihoodDelegate.BranchNodeOperation> branchNodeOperations,
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
        public final void simulate(final List<BranchNodeOperation> branchNodeOperations, final int rootNodeNumber) {

            setupStatistics();

            simulateRoot(rootNodeNumber);

            for (BranchNodeOperation operation : branchNodeOperations) {
                simulateNode(operation);
            }
        }

        @Override
        public final TreeTrait[] getTreeTraits() {
            return treeTraitHelper.getTreeTraits();
        }

        @Override
        public final TreeTrait getTreeTrait(String key) {
            return treeTraitHelper.getTreeTrait(key);
        }

        protected abstract void setupStatistics();

        protected abstract void simulateRoot(final int rootNumber);

        protected abstract void simulateNode(final BranchNodeOperation operation);

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
            this.diffusionModel = diffusionModel;

            sample = new double[dimNode * tree.getNodeCount()];

            tmpEpsilon = new double[dimTrait];
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
            System.arraycopy(sample, 0, trait, 0, length);

            return trait;
        }

        private double[] getTraitForNode(final NodeRef node) {

            assert simulationProcess != null;

            simulationProcess.cacheSimulatedTraits(null);

            double[] trait = new double[dimNode];
            System.arraycopy(sample, node.getNumber() * dimNode, trait, 0, dimNode);

            return trait;
        }

        private static double[][] getCholeskyOfVariance(double[][] precision) {
            Matrix variance = new SymmetricMatrix(precision).inverse();
            final double[][] cholesky;
            try {
                cholesky = new CholeskyDecomposition(variance).getL();
            } catch (IllegalDimension illegalDimension) {
                throw new RuntimeException("Attempted Cholesky decomposition on non-square matrix");
            }
            return cholesky;
        }

        @Override
        protected void setupStatistics() {
            double[][] diffusionPrecision = diffusionModel.getPrecisionmatrix();
            cholesky = getCholeskyOfVariance(diffusionPrecision); // TODO Cache
        }

        protected final double[] sample;

        protected final MultivariateDiffusionModel diffusionModel;

        protected double[][] cholesky;
        protected final double[] tmpEpsilon;
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

            this.likelihoodDelegate = likelihoodDelegate;
//            integrator = likelihoodDelegate.getIntegrator();

            final int partialLength = dataModel.getTipPartial(0).length; // TODO Need to generalize
            partialNodeBuffer = new double[partialLength];

//            tmpEpsilon = new double[dimTrait];
            tmpMean = new double[dimTrait];
        }

        @Override
        protected void simulateRoot(final int nodeIndex) {
            likelihoodDelegate.getPartial(nodeIndex, partialNodeBuffer);

            int offsetPartial = 0;
            int offsetSample = dimNode * nodeIndex;
            for (int trait = 0; trait < numTraits; ++trait) {

                final double nodePrecision = partialNodeBuffer[offsetPartial + dimTrait];

                if (Double.isInfinite(nodePrecision)) {
                    System.arraycopy(partialNodeBuffer, offsetPartial, sample, offsetSample, dimTrait);
                } else {

                    final double sqrtScale = Math.sqrt(1.0 / nodePrecision);

                    MultivariateNormalDistribution.nextMultivariateNormalCholesky(
                            partialNodeBuffer, offsetPartial, // input mean
                            cholesky, sqrtScale, // input variance
                            sample, offsetSample, // output sample
                            tmpEpsilon);
                }

                offsetSample += dimTrait;
                offsetPartial += (dimTrait + 1); // TODO Need to generalize;
            }
        }

        @Override
        protected void simulateNode(final BranchNodeOperation operation) {
            final int nodeIndex = operation.getNodeNumber();
            likelihoodDelegate.getPartial(nodeIndex, partialNodeBuffer);

            int offsetPartial = 0;
            int offsetSample = dimNode * nodeIndex;
            int offsetParent = dimNode * operation.getParentNumber();

            final double branchPrecision = 1.0 / operation.getBranchLength();
            
            for (int trait = 0; trait < numTraits; ++trait) {

                final double nodePrecision = partialNodeBuffer[offsetPartial + dimTrait];

                if (Double.isInfinite(nodePrecision)) {
                    System.arraycopy(partialNodeBuffer, offsetPartial, sample, offsetSample, dimTrait);
                } else {

                    final double totalPrecision = nodePrecision + branchPrecision;

                     for (int i = 0; i < dimTrait; ++i) {
                        tmpMean[i] = (nodePrecision * partialNodeBuffer[offsetPartial + i]
                                + branchPrecision * sample[offsetParent + i]) / totalPrecision;
                    }

                    final double sqrtScale = Math.sqrt(1.0 / totalPrecision);

                    MultivariateNormalDistribution.nextMultivariateNormalCholesky(
                            tmpMean, 0, // input mean
                            cholesky, sqrtScale, // input variance
                            sample, offsetSample, // output sample
                            tmpEpsilon);
                }

                offsetSample += dimTrait;
                offsetParent += dimTrait;
                offsetPartial += (dimTrait + 1); // TODO Need to generalize
            }
        }

//        private final MultivariateDiffusionModel diffusionModel;
//        private final ContinuousDiffusionIntegrator integrator;
        private final ContinuousDataLikelihoodDelegate likelihoodDelegate;
        private final double[] partialNodeBuffer;
//        private final double[] tmpEpsilon;
        private final double[] tmpMean;

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

            this.rootPrior = rootPrior;
        }

        @Override
        protected void simulateRoot(final int nodeIndex) {

            final double[] rootMean = rootPrior.getMean();
            final double sqrtScale = Math.sqrt(1.0 / rootPrior.getPseudoObservations());

            int offsetSample = dimNode * nodeIndex;
            for (int trait = 0; trait < numTraits; ++trait) {
                MultivariateNormalDistribution.nextMultivariateNormalCholesky(
                        rootMean, 0, // input meant
                        cholesky, sqrtScale,
                        sample, offsetSample,
                        tmpEpsilon
                );

                offsetSample += dimTrait;
            }
        }

        @Override
        protected void simulateNode(final BranchNodeOperation operation) {
            final int nodeIndex = operation.getNodeNumber();
            int offsetSample = dimNode * nodeIndex;
            int offsetParent = dimNode * operation.getParentNumber();

            final double branchLength =  operation.getBranchLength();

            if (branchLength == 0.0) {
                System.arraycopy(sample, offsetParent, sample, offsetSample, dimTrait * numTraits);
            } else {

                final double sqrtScale = Math.sqrt(branchLength);
                for (int trait = 0; trait < numTraits; ++trait) {
                    MultivariateNormalDistribution.nextMultivariateNormalCholesky(
                            sample, offsetParent,
                            cholesky, sqrtScale,
                            sample, offsetSample,
                            tmpEpsilon
                    );

                    offsetParent += dimTrait;
                    offsetSample += dimTrait;
                }
            }
        }

        private final ConjugateRootTraitPrior rootPrior;
    }

}
