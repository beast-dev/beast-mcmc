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

package dr.evomodel.treedatalikelihood.preorder;

import dr.evolution.tree.*;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.treedatalikelihood.ProcessOnTreeDelegate;
import dr.evomodel.treedatalikelihood.ProcessSimulation;
import dr.evomodel.treedatalikelihood.TreeTraversal;
import dr.evomodel.treedatalikelihood.continuous.*;
import dr.evomodel.treedatalikelihood.continuous.cdi.ContinuousDiffusionIntegrator;
import dr.inference.model.Model;
import dr.inference.model.ModelListener;
import dr.math.matrixAlgebra.*;
import org.ejml.data.DenseMatrix64F;

import java.util.List;
import java.util.Map;

import static dr.math.matrixAlgebra.missingData.MissingOps.*;

/**
 * ProcessSimulationDelegate - interface for a plugin delegate for data simulation on a tree.
 *
 * @author Andrew Rambaut
 * @author Marc Suchard
 * @version $Id$
 */
public interface ProcessSimulationDelegate extends ProcessOnTreeDelegate, TreeTraitProvider, ModelListener {

    void simulate(int[] operations, int operationCount, int rootNodeNumber);

    void setCallback(ProcessSimulation simulationProcess);

    int vectorizeNodeOperations(List<ProcessOnTreeDelegate.NodeOperation> nodeOperations, int[] operations);

    abstract class AbstractDelegate implements ProcessSimulationDelegate {

        AbstractDelegate(String name, Tree tree) {
            this.name = name;
            this.tree = tree;
            this.baseTree = getBaseTree(tree);
            constructTraits(treeTraitHelper);
        }

        protected abstract void constructTraits(Helper treeTraitHelper);

        @Override
        public final TreeTraversal.TraversalType getOptimalTraversalType() {
            return TreeTraversal.TraversalType.PRE_ORDER;
        }

        @Override
        public final void setCallback(ProcessSimulation simulationProcess) {
            this.simulationProcess = simulationProcess;
        }

        @Override
        public void simulate(final int[] operations, final int operationCount,
                             final int rootNodeNumber) {

            setupStatistics();

            simulateRoot(rootNodeNumber);

            int k = 0;
            for (int i = 0; i < operationCount; ++i) {
                simulateNode(
                        operations[k    ],
                        operations[k + 1],
                        operations[k + 2],
                        operations[k + 3],
                        operations[k + 4]
                );

                k += ContinuousDiffusionIntegrator.OPERATION_TUPLE_SIZE;
            }
        }

        private static Tree getBaseTree(Tree derived) {
            while (derived instanceof TransformableTree) {
                derived = ((TransformableTree) derived).getOriginalTree();
            }
            return derived;
        }

        static NodeRef getBaseNode(Tree derived, NodeRef node) {
            while (derived instanceof TransformableTree) {
                derived = ((TransformableTree) derived).getOriginalTree();
                node = ((TransformableTree) derived).getOriginalNode(node);
            }
            return node;
        }

        protected double getNormalization() {
            return 1.0;
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

        protected abstract void simulateNode(final int v0,
                                             final int v1,
                                             final int v2,
                                             final int v3,
                                             final int v4);

        final TreeTraitProvider.Helper treeTraitHelper = new Helper();

        ProcessSimulation simulationProcess = null;

        final Tree tree;
        final Tree baseTree;
        final String name;
    }

    abstract class AbstractContinuousTraitDelegate extends AbstractDelegate {

        final int dimTrait;
        final int numTraits;
        final int dimNode;

        final MultivariateDiffusionModel diffusionModel;
        final ContinuousTraitPartialsProvider dataModel;
        final ConjugateRootTraitPrior rootPrior;
        final RootProcessDelegate rootProcessDelegate;
        final ContinuousDataLikelihoodDelegate likelihoodDelegate;

        double[] diffusionVariance;
        DenseMatrix64F Vd;
        DenseMatrix64F Pd;

        double[][] cholesky;
        Map<PartiallyMissingInformation.HashedIntArray,
                ConditionalVarianceAndTransform> conditionalMap;

        AbstractContinuousTraitDelegate(String name,
                                        Tree tree,
                                        MultivariateDiffusionModel diffusionModel,
                                        ContinuousTraitPartialsProvider dataModel,
                                        ConjugateRootTraitPrior rootPrior,
                                        ContinuousRateTransformation rateTransformation,
                                        ContinuousDataLikelihoodDelegate likelihoodDelegate) {
            super(name, tree);

            dimTrait = likelihoodDelegate.getTraitDim();
            numTraits = likelihoodDelegate.getTraitCount();
            dimNode = dimTrait * numTraits;
            this.diffusionModel = diffusionModel;
            this.dataModel = dataModel;
            this.rateTransformation = rateTransformation;
            this.rootPrior = rootPrior;
            this.rootProcessDelegate = likelihoodDelegate.getRootProcessDelegate();
            this.likelihoodDelegate = likelihoodDelegate;

            diffusionModel.addModelListener(this);
        }

        @Override
        protected final double getNormalization() {
            return rateTransformation.getNormalization();
        }

        final private ContinuousRateTransformation rateTransformation;

        protected boolean isLoggable() {
            return true;
        }

        @Override
        public void modelChangedEvent(Model model, Object object, int index) {
            if (model == diffusionModel) {
                clearCache();
            } else {
                throw new IllegalArgumentException("Unknown model");
            }
        }

        @Override
        public void modelRestored(Model model) {

        }

        @Override
        protected void setupStatistics() {
            if (diffusionVariance == null) {
                double[][] diffusionPrecision = diffusionModel.getPrecisionmatrix();
                diffusionVariance = getVectorizedVarianceFromPrecision(diffusionPrecision);
                Vd = wrap(diffusionVariance, 0, dimTrait, dimTrait);
                Pd = new DenseMatrix64F(diffusionPrecision);
            }
            if (cholesky == null) {
                cholesky = getCholeskyOfVariance(diffusionVariance, dimTrait);
            }
        }

        void clearCache() {
            diffusionVariance = null;
            Vd = null;
            Pd = null;
            cholesky = null;
            conditionalMap = null;
        }

        static double[][] getCholeskyOfVariance(Matrix variance) {
            final double[][] cholesky;
            try {
                cholesky = new CholeskyDecomposition(variance).getL();
            } catch (IllegalDimension illegalDimension) {
                throw new RuntimeException("Attempted Cholesky decomposition on non-square matrix");
            }
            return cholesky;
        }

        static double[][] getCholeskyOfVariance(double[] variance, final int dim) {
            return CholeskyDecomposition.execute(variance, 0, dim);
        }

        private static double[] getVectorizedVarianceFromPrecision(double[][] precision) {
            return new SymmetricMatrix(precision).inverse().toVectorizedComponents();
        }

    }
}
