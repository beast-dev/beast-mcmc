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
import dr.evomodel.treedatalikelihood.continuous.*;
import dr.evomodel.treedatalikelihood.continuous.cdi.ContinuousDiffusionIntegrator;
import dr.evomodel.treedatalikelihood.continuous.cdi.PrecisionType;
import dr.inference.model.Model;
import dr.inference.model.ModelListener;
import dr.inference.model.Parameter;
import dr.math.MathUtils;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.math.matrixAlgebra.*;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static dr.evomodel.treedatalikelihood.continuous.cdi.ContinuousDiffusionIntegrator.Multivariate.*;

/**
 * ProcessSimulationDelegate - interface for a plugin delegate for data simulation on a tree.
 *
 * @author Andrew Rambaut
 * @author Marc Suchard
 * @version $Id$
 */
public interface ProcessSimulationDelegate extends ProcessOnTreeDelegate, TreeTraitProvider, ModelListener {

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

    abstract class AbstractDiscreteTraitDelegate extends AbstractDelegate {

        AbstractDiscreteTraitDelegate(String name, Tree tree) {
            super(name, tree);
        }
    }

    abstract class AbstractContinuousTraitDelegate extends AbstractDelegate {

        protected final int dimTrait;
        protected final int numTraits;
        protected final int dimNode;

        protected final MultivariateDiffusionModel diffusionModel;
        protected final ContinuousTraitDataModel dataModel;
        protected double[] diffusionVariance;
        protected DenseMatrix64F Vd;
        protected double[][] cholesky;
        protected Map<PartiallyMissingInformation.HashedIntArray,
                ConditionalOnPartiallyMissingTipsDelegate.ConditionalVarianceAndTranform> conditionalMap;

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
            this.dataModel = dataModel;

            diffusionModel.addModelListener(this);
        }

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
//                System.err.println("PDS.sS diffusionVariance");
                double[][] diffusionPrecision = diffusionModel.getPrecisionmatrix();
                diffusionVariance = getVectorizedVarianceFromPrecision(diffusionPrecision);
                Vd = wrap(diffusionVariance, 0, dimTrait, dimTrait);
            }
            if (cholesky == null) {
//                System.err.println("PDS.sS cholesky");
                cholesky = getCholeskyOfVariance(diffusionVariance, dimTrait); // TODO Cache
            }
        }

        public void clearCache() {
            diffusionVariance = null;
            Vd = null;
            cholesky = null;
            conditionalMap = null;
        }

        protected static double[][] getCholeskyOfVariance(Matrix variance) {
            final double[][] cholesky;
            try {
                cholesky = new CholeskyDecomposition(variance).getL();
            } catch (IllegalDimension illegalDimension) {
                throw new RuntimeException("Attempted Cholesky decomposition on non-square matrix");
            }
            return cholesky;
        }

        protected static double[][] getCholeskyOfVariance(double[] variance, final int dim) {
            return CholeskyDecomposition.execute(variance, 0, dim);
        }

        private static Matrix getVarianceFromPrecision(double[][] precision) {
            return new SymmetricMatrix(precision).inverse();
        }

        private static double[] getVectorizedVarianceFromPrecision(double[][] precision) {
            return new SymmetricMatrix(precision).inverse().toVectorizedComponents();
        }
    }

    class TipFullConditionalDistributionDelegate extends AbstractContinuousTraitDelegate {

        TipFullConditionalDistributionDelegate(String name, MultivariateTraitTree tree,
                                               MultivariateDiffusionModel diffusionModel,
                                               ContinuousTraitDataModel dataModel,
                                               ConjugateRootTraitPrior rootPrior,
                                               ContinuousRateTransformation rateTransformation,
                                               BranchRateModel rateModel,
                                               ContinuousDataLikelihoodDelegate likelihoodDelegate) {

            super(name, tree, diffusionModel, dataModel, rootPrior, rateTransformation, rateModel, likelihoodDelegate);
        }

        void constructTraits(Helper treeTraitHelper) {

            TreeTrait<MeanAndVariance> baseTrait = new TreeTrait<MeanAndVariance>() {

                public String getTraitName() { return "fcd." + name; }

                public Intent getIntent() { return Intent.NODE; }

                public Class getTraitClass() { return MeanAndVariance.class; }

                public MeanAndVariance getTrait(Tree t, NodeRef node) {
                    assert (tree == t);

                    return getTraitForNode(node);
                }

                public String getTraitString(Tree tree, NodeRef node) {
                    return getTrait(tree, node).toString();
                }

                public boolean getLoggable() { return isLoggable(); }
            };

            treeTraitHelper.addTrait(baseTrait);
        }

        private MeanAndVariance getTraitForNode(NodeRef node) {

            assert simulationProcess != null;
            assert node != null;

            if (nodeForLastCall != node) {
                // Re-simulate if calling for new node
                simulationProcess.modelChangedEvent(null, node, -1);
                nodeForLastCall = node;
            }

            simulationProcess.cacheSimulatedTraits(node);

            return cachedMeanAndVariance;
        }

        @Override
        protected void simulateRoot(int rootNumber) {
            System.err.println("computeRoot");
        }

        @Override
        protected void simulateNode(BranchNodeOperation operation) {
            System.err.println("computeNodes");
            cachedMeanAndVariance = new MeanAndVariance();
        }

        private NodeRef nodeForLastCall = null;
        private MeanAndVariance cachedMeanAndVariance;
    }

    class MeanAndVariance {
        double[] mean;
        Matrix variance;
    }

    abstract class AbstractRealizedContinuousTraitDelegate extends AbstractContinuousTraitDelegate {

        AbstractRealizedContinuousTraitDelegate(String name,
                                                MultivariateTraitTree tree,
                                                MultivariateDiffusionModel diffusionModel,
                                                ContinuousTraitDataModel dataModel,
                                                ConjugateRootTraitPrior rootPrior,
                                                ContinuousRateTransformation rateTransformation,
                                                BranchRateModel rateModel,
                                                ContinuousDataLikelihoodDelegate likelihoodDelegate) {
            super(name, tree, diffusionModel, dataModel, rootPrior, rateTransformation, rateModel, likelihoodDelegate);

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

            TreeTrait.DA tipTrait = new TreeTrait.DA() {

                public String getTraitName() { return getTipTraitName(name); }

                public Intent getIntent() { return Intent.WHOLE_TREE; }

                public double[] getTrait(Tree t, NodeRef node) {

                    assert t == tree;
                    return getTraitForAllTips();
                }
            };

            treeTraitHelper.addTrait(tipTrait);

            TreeTrait.DA tipPrecision = new TreeTrait.DA() {

                public String getTraitName() { return getTipPrecisionName(name); }

                public Intent getIntent() { return Intent.WHOLE_TREE; }

                public double[] getTrait(Tree t, NodeRef node) {

                    assert t == tree;
                    return getPrecisionForAllTips();
                }
            };

            treeTraitHelper.addTrait(tipPrecision);

        }

        public static String getTipTraitName(String name) {
            return "tip." + name;
        }

        public static String getTipPrecisionName(String name) {
            return "precision." + name;
        }

        private double[] getTraitForAllTips() {

            assert simulationProcess != null;

            simulationProcess.cacheSimulatedTraits(null);

            final int length = dimNode * tree.getExternalNodeCount();
            double[] trait = new double[length];
            System.arraycopy(sample, 0, trait, 0, length);

            return trait;
        }

        private double[] getPrecisionForAllTips() {

            assert simulationProcess != null;

            simulationProcess.cacheSimulatedTraits(null);

            final int length = tree.getExternalNodeCount();
            double[] precision = new double[length];

            Arrays.fill(precision, Double.POSITIVE_INFINITY); // TODO

            return precision;
        }

        private double[] getTraitForNode(final NodeRef node) {

            assert simulationProcess != null;

            simulationProcess.cacheSimulatedTraits(null);

            if (node == null) {
                return getTraitForAllTips();
            } else {
                double[] trait = new double[dimNode];
                System.arraycopy(sample, node.getNumber() * dimNode, trait, 0, dimNode);

                return trait;
            }

        }

        protected final double[] sample;
        protected final double[] tmpEpsilon;
    }

    class ConditionalOnTipsRealizedDelegate extends AbstractRealizedContinuousTraitDelegate {

        final protected int dimPartial;

        public ConditionalOnTipsRealizedDelegate(String name,
                                         MultivariateTraitTree tree,
                                         MultivariateDiffusionModel diffusionModel,
                                         ContinuousTraitDataModel dataModel,
                                         ConjugateRootTraitPrior rootPrior,
                                         ContinuousRateTransformation rateTransformation,
                                         BranchRateModel rateModel,
                                         ContinuousDataLikelihoodDelegate likelihoodDelegate) {
            super(name, tree, diffusionModel, dataModel, rootPrior, rateTransformation, rateModel, likelihoodDelegate);

            this.likelihoodDelegate = likelihoodDelegate;
            this.dimPartial = dimTrait + likelihoodDelegate.getPrecisionType().getMatrixLength(dimTrait);
            partialNodeBuffer = new double[numTraits * dimPartial];

            tmpMean = new double[dimTrait];
        }

        @Override
        protected void simulateRoot(final int nodeIndex) {
            likelihoodDelegate.getPartial(nodeIndex, partialNodeBuffer);

            int offsetPartial = 0;
            int offsetSample = dimNode * nodeIndex;
            for (int trait = 0; trait < numTraits; ++trait) {

                simulateTraitForRoot(offsetSample, offsetPartial);

                offsetSample += dimTrait;
                offsetPartial += dimPartial;
            }
        }

        protected void simulateTraitForRoot(final int offsetSample, final int offsetPartial) {

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
        }

        @Override
        protected void simulateNode(final BranchNodeOperation operation) {
            final int nodeIndex = operation.getNodeNumber();
            likelihoodDelegate.getPartial(nodeIndex, partialNodeBuffer);

            int offsetPartial = 0;
            int offsetSample = dimNode * nodeIndex;
            int offsetParent = dimNode * operation.getParentNumber();

            final double branchPrecision = 1.0 / operation.getBranchLength();

            final boolean isExternal = nodeIndex < tree.getExternalNodeCount();
            
            for (int trait = 0; trait < numTraits; ++trait) {

                simulateTraitForNode(offsetSample, offsetParent, offsetPartial, branchPrecision, isExternal);

                offsetSample += dimTrait;
                offsetParent += dimTrait;
                offsetPartial += dimPartial;
            }
        }

        protected void simulateTraitForNode(final int offsetSample,
                                             final int offsetParent,
                                             final int offsetPartial,
                                             final double branchPrecision,
                                            final boolean isExternal) {

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
        }

        protected final ContinuousDataLikelihoodDelegate likelihoodDelegate;
        protected final double[] partialNodeBuffer;
        protected final double[] tmpMean;
    }


    class MultivariateConditionalOnTipsRealizedDelegate extends ConditionalOnTipsRealizedDelegate {

        private static final boolean DEBUG = true;

        public MultivariateConditionalOnTipsRealizedDelegate(String name, MultivariateTraitTree tree,
                                                             MultivariateDiffusionModel diffusionModel,
                                                             ContinuousTraitDataModel dataModel,
                                                             ConjugateRootTraitPrior rootPrior,
                                                             ContinuousRateTransformation rateTransformation,
                                                             BranchRateModel rateModel,
                                                             ContinuousDataLikelihoodDelegate likelihoodDelegate) {
            super(name, tree, diffusionModel, dataModel, rootPrior, rateTransformation, rateModel, likelihoodDelegate);
        }

        @Override
        protected void simulateTraitForRoot(final int offsetSample, final int offsetPartial) {

            // TODO Bad programming -- should not need to know about internal layout
            // Layout, offset, dim
            // trait, 0, dT
            // precision, dT, dT * dT
            // variance, dT + dT * dT, dT * dT
            // scalar, dT + 2 * dT * dT, 1

            final double[][] cholesky = CholeskyDecomposition.execute(
                    partialNodeBuffer,
                    offsetPartial + dimTrait + dimTrait * dimTrait, // variance
                    dimTrait);  // TODO Reuse output memory

            MultivariateNormalDistribution.nextMultivariateNormalCholesky(
                    partialNodeBuffer, offsetPartial, // input mean
                    cholesky, 1.0, // input variance
                    sample, offsetSample, // output sample
                    tmpEpsilon);

            if (DEBUG) {
                System.err.println("Attempt to simulate root");
                final DenseMatrix64F mean = wrap(partialNodeBuffer, offsetPartial, dimTrait, 1);
                final DenseMatrix64F samp = wrap(sample, offsetSample, dimTrait, 1);
                final DenseMatrix64F V = wrap(partialNodeBuffer, offsetPartial + dimTrait + dimTrait * dimTrait, dimTrait, dimTrait);

                System.err.println("mean: " + mean);
                System.err.println("V: " + V);
                System.err.println("Ch:\n" + new Matrix(cholesky));
                System.err.println("sample: " + samp);
            }
        }

        @Override
        protected void simulateTraitForNode(final int offsetSample,
                                             final int offsetParent,
                                             final int offsetPartial,
                                             final double branchPrecision,
                                            final boolean isExternal) {

            if (isExternal) {
                simulateTraitForExternalNode(offsetSample, offsetParent, offsetPartial, branchPrecision);
            } else {
                simulateTraitForInternalNode(offsetSample, offsetParent, offsetPartial, branchPrecision);
            }

        }


        boolean debugBreak = false;

        private void simulateTraitForExternalNode(final int offsetSample,
                                             final int offsetParent,
                                             final int offsetPartial,
                                             final double branchPrecision) {

            final DenseMatrix64F P0 = wrap(partialNodeBuffer, offsetPartial + dimTrait, dimTrait, dimTrait);

            if (DEBUG) {
                System.err.println("sTFEN");
                System.err.println("P0: " + P0);
            }

            final int missingCount = countFiniteDiagonals(P0);

            if (missingCount == 0) {
                System.arraycopy(partialNodeBuffer, offsetPartial, sample, offsetSample, dimTrait);
            } else {

                if (DEBUG) {
                    System.err.println("Missing count: " + missingCount);
                    debugBreak = true;
                }
                simulateTraitForInternalNode(offsetSample, offsetParent, offsetPartial, branchPrecision);

//                for (int i = 0; i < missingCount; ++i) {
//                    MathUtils.nextGaussian();
//                }
            }

//            System.exit(-1);
        }

        private void simulateTraitForInternalNode(final int offsetSample,
                                             final int offsetParent,
                                             final int offsetPartial,
                                             final double branchPrecision) {

            final WrappedVector M0 = new WrappedVector(partialNodeBuffer, offsetPartial, dimTrait);
            final DenseMatrix64F P0 = wrap(partialNodeBuffer, offsetPartial + dimTrait, dimTrait, dimTrait);
//            final DenseMatrix64F V0 = wrap(partialNodeBuffer, offsetPartial + dimTrait + dimTrait * dimTrait, dimTrait, dimTrait);

            final WrappedVector M1 = new WrappedVector(sample, offsetParent, dimTrait);
            final DenseMatrix64F P1 = new DenseMatrix64F(dimTrait, dimTrait);
            final DenseMatrix64F V1 = new DenseMatrix64F(dimTrait, dimTrait);

            CommonOps.scale(1.0 / branchPrecision, Vd, V1);
            final ContinuousDiffusionIntegrator.Multivariate.InversionResult c1 = safeInvert(V1, P1, false);

            final WrappedVector M2 = new WrappedVector(tmpMean, 0, dimTrait);
            final DenseMatrix64F P2 = new DenseMatrix64F(dimTrait, dimTrait);
            final DenseMatrix64F V2 = new DenseMatrix64F(dimTrait, dimTrait);

            CommonOps.add(P0, P1, P2);
            final ContinuousDiffusionIntegrator.Multivariate.InversionResult c2 = safeInvert(P2, V2, false);
            weightedAverage(M0, P0, M1, P1, M2, V2, dimTrait);

            double[][] C2 = getCholeskyOfVariance(V2.getData(), dimTrait);

            MultivariateNormalDistribution.nextMultivariateNormalCholesky(
                    M2.getBuffer(), 0, // input mean
                    C2, 1.0, // input variance
                    sample, offsetSample, // output sample
                    tmpEpsilon);

            if (DEBUG) {
                System.err.println("sTFIN");
                System.err.println("M0: " + M0);
                System.err.println("P0: " + P0);
                System.err.println("M1: " + M1);
                System.err.println("P1: " + P1);
                System.err.println("M2: " + M2);
                System.err.println("V2: " + V2);
                System.err.println("C2: " + new Matrix(C2));
                System.err.println("SS: " + new WrappedVector(sample, offsetSample, dimTrait));
                System.err.println("");

                if (debugBreak) {
                    debugBreak = false;
//                    System.exit(-1);
                }
//                System.exit(-1);
            }

        }
    }

    class ConditionalOnPartiallyMissingTipsDelegate extends ConditionalOnTipsRealizedDelegate {

        @Deprecated
        public ConditionalOnPartiallyMissingTipsDelegate(String name, MultivariateTraitTree tree,
                                                         MultivariateDiffusionModel diffusionModel,
                                                         ContinuousTraitDataModel dataModel,
                                                         ConjugateRootTraitPrior rootPrior,
                                                         ContinuousRateTransformation rateTransformation,
                                                         BranchRateModel rateModel,
                                                         ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                                         Parameter missingParameter) {

            this(name, tree, diffusionModel, dataModel, rootPrior, rateTransformation, rateModel, likelihoodDelegate,
                    new PartiallyMissingInformation(tree.getExternalNodeCount(), dataModel.getTraitCount(),
                            dataModel.getTraitDimension(), missingParameter));
        }

        public ConditionalOnPartiallyMissingTipsDelegate(String name, MultivariateTraitTree tree,
                                                         MultivariateDiffusionModel diffusionModel,
                                                         ContinuousTraitDataModel dataModel,
                                                         ConjugateRootTraitPrior rootPrior,
                                                         ContinuousRateTransformation rateTransformation,
                                                         BranchRateModel rateModel,
                                                         ContinuousDataLikelihoodDelegate likelihoodDelegate) {

            this(name, tree, diffusionModel, dataModel, rootPrior, rateTransformation, rateModel, likelihoodDelegate,
                    new PartiallyMissingInformation(tree, dataModel, likelihoodDelegate));
        }

        public ConditionalOnPartiallyMissingTipsDelegate(String name, MultivariateTraitTree tree,
                                                         MultivariateDiffusionModel diffusionModel,
                                                         ContinuousTraitDataModel dataModel,
                                                         ConjugateRootTraitPrior rootPrior,
                                                         ContinuousRateTransformation rateTransformation,
                                                         BranchRateModel rateModel,
                                                         ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                                         PartiallyMissingInformation missingInformation) {

            super(name, tree, diffusionModel, dataModel, rootPrior, rateTransformation, rateModel, likelihoodDelegate);
            this.missingInformation = missingInformation;

            assert (dataModel.getPrecisionType() == PrecisionType.FULL);
        }

        @Override
        protected boolean isLoggable() {
            return false;
        }

        final private PartiallyMissingInformation missingInformation;

        @Override
        protected void simulateNode(final BranchNodeOperation operation) {
            final int nodeIndex = operation.getNodeNumber();
            likelihoodDelegate.getPartial(nodeIndex, partialNodeBuffer);

            int offsetPartial = 0;
            int offsetSample = dimNode * nodeIndex;
            int offsetParent = dimNode * operation.getParentNumber();

            final boolean isExternal = nodeIndex < tree.getExternalNodeCount();

            final double branchPrecision = 1.0 / operation.getBranchLength();

            for (int trait = 0; trait < numTraits; ++trait) {

                final double nodePrecision = partialNodeBuffer[offsetPartial + dimTrait];  // TODO PrecisionType.FULL

                if (!isExternal) {

                    simulateTraitForNode(offsetSample, offsetParent, offsetPartial, nodePrecision, isExternal);

                } else { // Is external

                    // Copy tip values into sample
                    System.arraycopy(partialNodeBuffer, offsetPartial, sample, offsetSample, dimTrait);

                    if (missingInformation.isPartiallyMissing(nodeIndex, trait)) {

                        PartiallyMissingInformation.HashedIntArray missingIndices =
                                missingInformation.getMissingIndices(nodeIndex, trait);

                        final int numMissing = missingIndices.getLength();
                        final int numNotMissing = missingIndices.getComplementLength();

                        assert (numMissing + numNotMissing == dimTrait);

                        ConditionalVarianceAndTranform transform;
                        try {
                            transform = conditionalMap.get(missingIndices);
                        } catch (NullPointerException nep) {
//                            System.err.println("Make CVT");
                            transform =
//                                    new ConditionalVarianceAndTranform(diffusionVariance,
//                                    missingIndices.getArray(),
//                                    missingIndices.getComplement());
                            null;

                            if (conditionalMap == null) {
                                conditionalMap = new HashMap<PartiallyMissingInformation.HashedIntArray,
                                        ConditionalVarianceAndTranform>();
                            }
                            conditionalMap.put(missingIndices, transform);
                        }
                        // TODO Must clear cache

//                        ConditionalVarianceAndTranform transform =
//                                new ConditionalVarianceAndTranform(diffusionVariance,
//                                        missingIndices.getArray(),
//                                        missingIndices.getComplement());

                        // TODO PrecisionType.FULL

                        final double[] conditionalMean = transform.getConditionalMean(sample, offsetSample,
                                sample, offsetParent);
                        final double[][] conditionalCholesky = transform.getConditionalCholesky();

                        final double sqrtScale = Math.sqrt(1.0 / branchPrecision);

                        MultivariateNormalDistribution.nextMultivariateNormalCholesky(
                                conditionalMean, 0, // input mean
                                conditionalCholesky, sqrtScale, // input variance
                                tmpMean, 0, // output sample
                                transform.getTemporageStorage());

                        for (int i = 0; i < numMissing; ++i) {
                            sample[offsetSample + missingIndices.get(i)] = tmpMean[i];
                        }

//                        System.err.println("mean:\n" + new Vector(conditionalMean));
//                        System.err.println("cholesky:\n" + new Matrix(conditionalCholesky));
//                        System.err.println("sS: " + sqrtScale);
//                        System.err.println("cMean\n" + new Vector(tmpMean));
//                        System.err.println("");
//                        System.err.println("");
                    }
                }

                offsetSample += dimTrait;
                offsetParent += dimTrait;
                offsetPartial += (dimTrait + 1); // TODO PrecisionType.FULL
            }
        }

        class ConditionalVarianceAndTranform {

            /**
             * For partially observed tips: (y_1, y_2)^t \sim N(\mu, \Sigma) where
             *
             *      \mu = (\mu_1, \mu_2)^t
             *      \Sigma = ((\Sigma_{11}, \Sigma_{12}), (\Sigma_{21}, \Sigma_{22})^t
             *
             * then  y_1 | y_2 \sim N (\bar{\mu}, \bar{\Sigma}), where
             *
             *      \bar{\mu} = \mu_1 + \Sigma_{12}\Sigma_{22}^{-1}(y_2 - \mu_2), and
             *      \bar{\Sigma} = \Sigma_{11} - \Sigma_{12}\Sigma_{22}^1\Sigma{21}
             *
             */

            final private double[][] cholesky;
            final private Matrix affineTransform;
            final int[] missingIndices;
            final int[] notMissingIndices;
            final double[] tempStorage;

            final int numMissing;
            final int numNotMissing;

            private static final boolean DEBUG = false;

            ConditionalVarianceAndTranform(final Matrix variance, final int[] missingIndices, final int[] notMissingIndices) {

                assert (missingIndices.length + notMissingIndices.length == variance.rows());
                assert (missingIndices.length + notMissingIndices.length == variance.columns());

                this.missingIndices = missingIndices;
                this.notMissingIndices = notMissingIndices;

                if (DEBUG) {
                    System.err.println("variance:\n" + variance);
                }

                Matrix S12S22Inv = null;
                Matrix Sbar = null;

                try {

                    Matrix S22 = variance.extractRowsAndColumns(notMissingIndices, notMissingIndices);
                    if (DEBUG) {
                        System.err.println("S22:\n" + S22);
                    }

                    Matrix S22Inv = S22.inverse();
                    if (DEBUG) {
                        System.err.println("S22Inv:\n" + S22Inv);
                    }

                    Matrix S12 = variance.extractRowsAndColumns(missingIndices, notMissingIndices);
                    if (DEBUG) {
                        System.err.println("S12:\n" + S12);
                    }

                    S12S22Inv = S12.product(S22Inv);
                    if (DEBUG) {
                        System.err.println("S12S22Inv:\n" + S12S22Inv);
                    }

                    Matrix S12S22InvS21 = S12S22Inv.productWithTransposed(S12);
                    if (DEBUG) {
                        System.err.println("S12S22InvS21:\n" + S12S22InvS21);
                    }

                    Sbar = variance.extractRowsAndColumns(missingIndices, missingIndices);
                    Sbar.decumulate(S12S22InvS21);
                    if (DEBUG) {
                        System.err.println("Sbar:\n" + Sbar);
                    }

                } catch (IllegalDimension illegalDimension) {
                    illegalDimension.printStackTrace();
                }

                this.affineTransform = S12S22Inv;
                this.cholesky = getCholeskyOfVariance(Sbar);
                this.tempStorage = new double[missingIndices.length];

                this.numMissing = missingIndices.length;
                this.numNotMissing = notMissingIndices.length;

            }

            double[] getConditionalMean(final double[] y, final int offsetY,
                                        final double[] mu, final int offsetMu) {

                double[] muBar = new double[numMissing];

                double[] shift = new double[numNotMissing];
                for (int i = 0; i < numNotMissing; ++i) {
                    final int noti = notMissingIndices[i];
                    shift[i] = y[offsetY + noti] - mu[offsetMu + noti];
                }

                for (int i = 0; i < numMissing; ++i) {
                    double delta = 0.0;
                    for (int k = 0; k < numNotMissing; ++k) {
                        delta += affineTransform.component(i, k) * shift[k];
                    }

                    muBar[i] = mu[offsetMu + missingIndices[i]] + delta;
                }

                return muBar;
            }

            void scatterResult(final double[] source, final int offsetSource,
                               final double[] destination, final int offsetDestination) {
                for (int i = 0; i < numMissing; ++i) {
                    destination[offsetDestination + missingIndices[i]] = source[offsetSource + i];
                }
            }

            double[][] getConditionalCholesky() {
                return cholesky;
            }

            Matrix getAffineTransform() {
                return affineTransform;
            }

            double[] getTemporageStorage() {
                return tempStorage;
            }
        }

//        @Override
//        protected void doDraw(final int nodeIndex, final double[] mean, final double totalPrecision, final int offsetSample) {
//
////            System.err.println(nodeIndex);
////            System.exit(-1);
//
//            final boolean isExternal = nodeIndex < tree.getExternalNodeCount();
//
//            if (isExternal) {
//
////                final boolean[] missing = dataModel.getPartiallyMissing(nodeIndex);
//
//                throw new RuntimeException("Not yet implemented");
//            } else {
//                super.doDraw(nodeIndex, mean, 0, totalPrecision, offsetSample);
//            }
//        }

        public static final String PARTIAL = "partial";

        public static String getPartiallyMissingTraitName(final String traitName) {
            return PARTIAL + "." + traitName;
        }
    }

    class UnconditionalOnTipsDelegate extends AbstractRealizedContinuousTraitDelegate {

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
                        rootMean, 0, // input mean
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

    class NotHere {
//        CholeskyDecompositionImpl impl = new CholeskyDecompositionCommon_D64()
//        impl.



        public static DenseMatrix64F wrap(final double[] source, final int offset,
                                          final int numRows, final int numCols,
                                          final double[] buffer) {
            System.arraycopy(source, offset, buffer, 0, numRows * numCols);
            return DenseMatrix64F.wrap(numRows, numCols, buffer);
        }

        public NotHere() {
            play();
        }

        private void play() {

            DenseMatrix64F A = DenseMatrix64F.wrap(2, 2, new double[4]);
//            DenseMatrix64F Ainv = DenseMatrix64F.wrap(2, 2, new double[4]);

//            LinearSolver<DenseMatrix64F> solver = LinearSolverFactory.qr(A.numRows,A.numCols);
//            if (!solver.setA(A)) {
//                throw new RuntimeException("Singular matrix");
//            }
//
//            solver.invert(Ainv);

            CommonOps.invert(A);
            System.err.println(A.toString());


//            if( !solver.setA(A) ) {
//                throw new IllegalArgument("Singular matrix");
//            }
//
//            if( solver.quality() <= 1e-8 )
//                throw new IllegalArgument("Nearly singular matrix");
//
//            solver.solve(b,x);
        }
    }

}
