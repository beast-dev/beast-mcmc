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
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Model;
import dr.inference.model.ModelListener;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.math.matrixAlgebra.*;
import dr.math.matrixAlgebra.missingData.InversionResult;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import java.util.Arrays;
import java.util.HashMap;
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

    void simulate(SimulationTreeTraversal treeTraversal,
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
        public void simulate(final SimulationTreeTraversal treeTraversal,
                                   final int rootNodeNumber) {

            final List<BranchNodeOperation> branchNodeOperations = treeTraversal.getBranchNodeOperations();
            final double normalization = getNormalization();

            setupStatistics();

            simulateRoot(rootNodeNumber);

            for (BranchNodeOperation operation : branchNodeOperations) {
                simulateNode(operation, normalization);
            }
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

        protected abstract void simulateNode(final BranchNodeOperation operation, final double branchNormalization);

        protected abstract void simulateNode(final NodeOperation operation);

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
        protected final ContinuousTraitPartialsProvider dataModel;
        protected final ConjugateRootTraitPrior rootPrior;
        protected final RootProcessDelegate rootProcessDelegate;

        protected double[] diffusionVariance;
        protected DenseMatrix64F Vd;
        protected DenseMatrix64F Pd;

        protected double[][] cholesky;
        protected Map<PartiallyMissingInformation.HashedIntArray,
                ConditionalOnPartiallyMissingTipsRealizedDelegate.ConditionalVarianceAndTranform> conditionalMap;

        AbstractContinuousTraitDelegate(String name,
                                        Tree tree,
                                        MultivariateDiffusionModel diffusionModel,
                                        ContinuousTraitPartialsProvider dataModel,
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
            this.rateTransformation = rateTransformation;
            this.rootPrior = rootPrior;
            this.rootProcessDelegate = likelihoodDelegate.getRootProcessDelegate();

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
//                System.err.println("PDS.sS diffusionVariance");
                double[][] diffusionPrecision = diffusionModel.getPrecisionmatrix();
                diffusionVariance = getVectorizedVarianceFromPrecision(diffusionPrecision);
                Vd = wrap(diffusionVariance, 0, dimTrait, dimTrait);
                Pd = new DenseMatrix64F(diffusionPrecision);
            }
            if (cholesky == null) {
//                System.err.println("PDS.sS cholesky");
                cholesky = getCholeskyOfVariance(diffusionVariance, dimTrait); // TODO Cache
            }
        }

        public void clearCache() {
            diffusionVariance = null;
            Vd = null;
            Pd = null;
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

    abstract class AbstractValuesViaFullConditionalDelegate extends TipFullConditionalDistributionDelegate {

        final private PartiallyMissingInformation missingInformation;

        protected boolean isLoggable() {
                    return false;
                }

        public AbstractValuesViaFullConditionalDelegate(String name, Tree tree,
                                                           MultivariateDiffusionModel diffusionModel,
                                                           ContinuousTraitPartialsProvider dataModel,
                                                           ConjugateRootTraitPrior rootPrior,
                                                           ContinuousRateTransformation rateTransformation,
                                                           BranchRateModel rateModel,
                                                           ContinuousDataLikelihoodDelegate likelihoodDelegate) {
            super(name, tree, diffusionModel, dataModel, rootPrior, rateTransformation, rateModel, likelihoodDelegate);
            missingInformation = new PartiallyMissingInformation(tree, dataModel, likelihoodDelegate);
        }
        
        protected double[] getTraitForNode(NodeRef node) {

            assert simulationProcess != null;
            assert node != null;

            final int nodeBuffer = likelihoodDelegate.getActiveNodeIndex(node.getNumber());

            if (node.getNumber() >= tree.getExternalNodeCount()) {   // Not external node
                return new double[0];
//                return new MeanAndVariance(new double[0]);
            }

            double[] conditionalNodeBuffer = null; //new double[dimPartial * numTraits];
            likelihoodDelegate.getPostOrderPartial(node.getNumber(), partialNodeBuffer);

            final double[] sample = new double[dimTrait * numTraits];

            int partialOffset = 0;
            int sampleOffset = 0;

            for (int trait = 0; trait < numTraits; ++trait) {
                if (missingInformation.isPartiallyMissing(node.getNumber(), trait)) {
                    if (conditionalNodeBuffer == null) {
                        conditionalNodeBuffer = new double[dimPartial * numTraits];

                        simulationProcess.cacheSimulatedTraits(node);
                        likelihoodDelegate.getPreOrderPartial(node.getNumber(), conditionalNodeBuffer);
                    }

                    System.err.println("Missing tip = " + node.getNumber() + " (" + nodeBuffer + "), trait = " + trait);

                    final WrappedVector preMean = new WrappedVector.Raw(conditionalNodeBuffer, partialOffset, dimTrait);
                    final DenseMatrix64F preVar = wrap(conditionalNodeBuffer, partialOffset + dimTrait + dimTrait * dimTrait, dimTrait, dimTrait);

                    final WrappedVector postObs = new WrappedVector.Raw(partialNodeBuffer, partialOffset, dimTrait);

                    System.err.println("post: " + postObs);
                    System.err.println("pre : " + preMean);
                    System.err.println("V: " + preVar);

                    if (missingInformation.isCompletelyMissing(node.getNumber(), trait)) {

                    } else {

                        final PartiallyMissingInformation.HashedIntArray intArray =
                                missingInformation.getMissingIndices(node.getNumber(), trait);
                        final int[] missing = intArray.getArray();
                        final int[] observed = intArray.getComplement();

                        ConditionalOnPartiallyMissingTipsRealizedDelegate.ConditionalVarianceAndTranform2 transform =
                                new ConditionalOnPartiallyMissingTipsRealizedDelegate.ConditionalVarianceAndTranform2(
                                        preVar, missing, observed
                                );

                        final WrappedVector cM = transform.getConditionalMean(
                                partialNodeBuffer, partialOffset,      // Tip value
                                conditionalNodeBuffer, partialOffset); // Mean value

                        computeValueWithMissing(cM, // input mean
                                transform.getConditionalCholesky(), // input variance,
                                new WrappedVector.Indexed(sample, sampleOffset, missing, missing.length), // output sample
                                transform.getTemporageStorage());

                        System.err.println("cM: " + cM);
                        System.err.println("CV: " + transform.getConditionalVariance());
                        System.err.println("value: " + new WrappedVector.Raw(sample, sampleOffset, dimTrait));
                    }

                } else {
                    computeValueWithNoMissing(partialNodeBuffer, partialOffset, sample, sampleOffset, dimTrait);
                }

                partialOffset += dimPartial;
                sampleOffset += dimTrait;
            }

            return sample;
//            return new MeanAndVariance(sample);
        }

        abstract protected void computeValueWithNoMissing(final double[] mean, final int meanOffset,
                                             final double[] output, final int outputOffset,
                                             final int dim);
        
        abstract protected void computeValueWithMissing(final WrappedVector mean,
                                               final double[][] cholesky,
                                               final WrappedVector output,
                                               final double[] buffer);
    }

    class TipRealizedValuesViaFullConditionalDelegate extends AbstractValuesViaFullConditionalDelegate {

        public TipRealizedValuesViaFullConditionalDelegate(String name, Tree tree,
                                                           MultivariateDiffusionModel diffusionModel,
                                                           ContinuousTraitPartialsProvider dataModel,
                                                           ConjugateRootTraitPrior rootPrior,
                                                           ContinuousRateTransformation rateTransformation,
                                                           BranchRateModel rateModel,
                                                           ContinuousDataLikelihoodDelegate likelihoodDelegate) {
            super(name, tree, diffusionModel, dataModel, rootPrior, rateTransformation, rateModel, likelihoodDelegate);
        }
        
        public static String getTraitName(String name) { return "tipSample." + name; }

        @Override
        protected void computeValueWithNoMissing(final double[] mean, final int meanOffset,
                                             final double[] output, final int outputOffset,
                                             final int dim) {
            System.arraycopy(mean, meanOffset, output, outputOffset, dim);
        }

        @Override
        protected void computeValueWithMissing(final WrappedVector mean,
                                               final double[][] cholesky,
                                               final WrappedVector output,
                                               final double[] buffer) {
            MultivariateNormalDistribution.nextMultivariateNormalCholesky(
                    mean, // input mean
                    cholesky, 1.0, // input variance
                    output, // output sample
                    buffer);
        }
    }

    class TipGradientViaFullConditionalDelegate extends TipFullConditionalDistributionDelegate {

//        final private PartiallyMissingInformation missingInformation;

         public TipGradientViaFullConditionalDelegate(String name, Tree tree,
                                                            MultivariateDiffusionModel diffusionModel,
                                                            ContinuousTraitDataModel dataModel,
                                                            ConjugateRootTraitPrior rootPrior,
                                                            ContinuousRateTransformation rateTransformation,
                                                            BranchRateModel rateModel,
                                                            ContinuousDataLikelihoodDelegate likelihoodDelegate) {
             super(name, tree, diffusionModel, dataModel, rootPrior, rateTransformation, rateModel, likelihoodDelegate);
//             missingInformation = new PartiallyMissingInformation(tree, dataModel, likelihoodDelegate);

             if (likelihoodDelegate.getPrecisionType() != PrecisionType.SCALAR) {
                 throw new RuntimeException("Tip gradients are not implemented for '" +
                         likelihoodDelegate.getPrecisionType().toString() + "' likelihoods");
             }
         }

         public static String getTraitName(String name) { return "grad." + name; }

         @Override
         protected double[] getTraitForNode(NodeRef node) {

             final double[] fullConditionalPartial = super.getTraitForNode(node);

             final double[] postOrderPartial = new double[dimPartial * numTraits];
             cdi.getPostOrderPartial(likelihoodDelegate.getActiveNodeIndex(node.getNumber()), postOrderPartial);

             final MatrixParameterInterface precision = diffusionModel.getPrecisionParameter();

             final double[] gradient = new double[dimTrait * numTraits];

             if (numTraits > 1) {
                 throw new RuntimeException("Not yet implemented");
             }

             final double scale = fullConditionalPartial[dimTrait];

             for (int i = 0; i < dimTrait; ++i) {

                 double sum = 0.0;
                 for (int j = 0; j < dimTrait; ++j) {
                     sum += (fullConditionalPartial[j] - postOrderPartial[j]) * scale *
                             precision.getParameterValue(i * dimTrait + j);
                 }

                 gradient[i] = sum;
             }

             return gradient;
         }
     }


    class TipFullConditionalDistributionDelegate extends AbstractContinuousTraitDelegate {

        public TipFullConditionalDistributionDelegate(String name, Tree tree,
                                               MultivariateDiffusionModel diffusionModel,
                                                      ContinuousTraitPartialsProvider dataModel,
                                               ConjugateRootTraitPrior rootPrior,
                                               ContinuousRateTransformation rateTransformation,
                                               BranchRateModel rateModel,
                                               ContinuousDataLikelihoodDelegate likelihoodDelegate) {

            super(name, tree, diffusionModel, dataModel, rootPrior, rateTransformation, rateModel, likelihoodDelegate);
            buffer = new MeanAndVariance[tree.getExternalNodeCount()];
            this.likelihoodDelegate = likelihoodDelegate;
            this.cdi = likelihoodDelegate.getIntegrator();

            this.dimPartial = dimTrait + likelihoodDelegate.getPrecisionType().getMatrixLength(dimTrait);
            partialNodeBuffer = new double[numTraits * dimPartial];

        }

        protected boolean isLoggable() {
                    return false;
                }

        protected final ContinuousDataLikelihoodDelegate likelihoodDelegate;
        protected final ContinuousDiffusionIntegrator cdi;
        protected final int dimPartial;
        protected final double[] partialNodeBuffer;

        public static String getTraitName(String name) { return "fcd." + name; }

        protected String delegateGetTraitName() {  return getTraitName(name); }

        protected Class delegateGetTraitClass() { return double[].class; }

        void constructTraits(Helper treeTraitHelper) {

            TreeTrait.DA baseTrait = new TreeTrait.DA() {

                public String getTraitName() { return delegateGetTraitName(); }

                public Intent getIntent() { return Intent.NODE; }

                public Class getTraitClass() { return delegateGetTraitClass(); }

                public double[] getTrait(Tree t, NodeRef node) {
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

        protected double[] getTraitForNode(NodeRef node) {
//        private double[] getTraitForNode(NodeRef node) {

            assert simulationProcess != null;
            assert node != null;

            simulationProcess.cacheSimulatedTraits(node);

            double[] partial = new double[dimPartial * numTraits];
            cdi.getPreOrderPartial(likelihoodDelegate.getActiveNodeIndex(node.getNumber()), partial);

            return partial;
//            return new MeanAndVariance(partial);
        }

        @Override
        protected void simulateRoot(int rootNumber) {

            if (DEBUG) {
                System.err.println("computeRoot");
            }

            final DenseMatrix64F diffusion = new DenseMatrix64F(likelihoodDelegate.getDiffusionModel().getPrecisionmatrix());

            // Copy post-order root prior to pre-order

            final double[] priorBuffer = partialNodeBuffer;
            final double[] rootBuffer = new double[priorBuffer.length];

            cdi.getPostOrderPartial(rootProcessDelegate.getPriorBufferIndex(), partialNodeBuffer); // No double-buffering

            int offset = 0;
            for (int trait = 0; trait < numTraits; ++trait) {
                // Copy mean
                System.arraycopy(priorBuffer, offset, rootBuffer, offset, dimTrait);

                final DenseMatrix64F Pp = wrap(priorBuffer, offset + dimTrait, dimTrait, dimTrait);
                final DenseMatrix64F Pr = new DenseMatrix64F(dimTrait, dimTrait);
                CommonOps.mult(Pp, diffusion, Pr);

                unwrap(Pr, rootBuffer, offset + dimTrait);

                offset += dimPartial;
            }

            cdi.setPreOrderPartial(likelihoodDelegate.getActiveNodeIndex(rootNumber), rootBuffer);

            if (DEBUG) {
                System.err.println("Root: " + new WrappedVector.Raw(rootBuffer, 0, rootBuffer.length));
                System.err.println("");
            }
        }

        @Override
        protected void simulateNode(BranchNodeOperation operation, double branchNormalization) {
            throw new RuntimeException("Not implemented");
        }

        @Override
        protected void simulateNode(NodeOperation operation) {

            cdi.updatePreOrderPartial(
                    likelihoodDelegate.getActiveNodeIndex(operation.getNodeNumber()),
                    likelihoodDelegate.getActiveNodeIndex(operation.getLeftChild()),
                    likelihoodDelegate.getActiveMatrixIndex(operation.getLeftChild()),
                    likelihoodDelegate.getActiveNodeIndex(operation.getRightChild()),
                    likelihoodDelegate.getActiveMatrixIndex(operation.getRightChild())
            );
        }

        @Override
        public final void simulate(final SimulationTreeTraversal treeTraversal,
                                   final int rootNodeNumber) {


            final List<NodeOperation> nodeOperations = treeTraversal.getNodeOperations();
            setupStatistics();

            simulateRoot(rootNodeNumber);

            for (NodeOperation operation : nodeOperations) {
                simulateNode(operation);
            }

            if (DEBUG) {
                System.err.println("END OF PRE-ORDER");
            }
        }

        final private MeanAndVariance[] buffer;

        private NodeRef nodeForLastCall = null;
        private MeanAndVariance cachedMeanAndVariance;

        private static final boolean DEBUG = false;

    }

    class MeanAndVariance {
        double[] mean;
        Matrix variance;

        public MeanAndVariance(double[] mean) {
            this.mean = mean;
        }

        public double[] getMean() { return mean; }

        public Matrix getVariance() { return variance; }
    }

    abstract class AbstractRealizedContinuousTraitDelegate extends AbstractContinuousTraitDelegate {

        AbstractRealizedContinuousTraitDelegate(String name,
                                                Tree tree,
                                                MultivariateDiffusionModel diffusionModel,
                                                ContinuousTraitPartialsProvider dataModel,
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

//                    assert t == tree; // Does not hold for transformed trees
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

        static final private boolean DEBUG = false;

        final protected int dimPartial;

        public ConditionalOnTipsRealizedDelegate(String name,
                                         Tree tree,
                                         MultivariateDiffusionModel diffusionModel,
                                         ContinuousTraitPartialsProvider dataModel,
                                         ConjugateRootTraitPrior rootPrior,
                                         ContinuousRateTransformation rateTransformation,
                                         BranchRateModel rateModel,
                                         ContinuousDataLikelihoodDelegate likelihoodDelegate) {
            super(name, tree, diffusionModel, dataModel, rootPrior, rateTransformation, rateModel, likelihoodDelegate);

            this.likelihoodDelegate = likelihoodDelegate;
            this.dimPartial = dimTrait + likelihoodDelegate.getPrecisionType().getMatrixLength(dimTrait);
            partialNodeBuffer = new double[numTraits * dimPartial];
            partialPriorBuffer = new double[numTraits * dimPartial];

            tmpMean = new double[dimTrait];
        }

        @Override
        protected void simulateRoot(final int nodeIndex) {

            likelihoodDelegate.getIntegrator().getPostOrderPartial(
                    likelihoodDelegate.getRootProcessDelegate().getPriorBufferIndex(),
                    partialPriorBuffer);

            likelihoodDelegate.getPostOrderPartial(nodeIndex, partialNodeBuffer);

            if (DEBUG) {
                System.err.println("Simulate root node " + nodeIndex);
            }

            int offsetPartial = 0;
            int offsetSample = dimNode * nodeIndex;
            for (int trait = 0; trait < numTraits; ++trait) {

                simulateTraitForRoot(offsetSample, offsetPartial);

                offsetSample += dimTrait;
                offsetPartial += dimPartial;
            }
        }

        protected void simulateTraitForRoot(final int offsetSample, final int offsetPartial) {

            final double rootPrec = partialNodeBuffer[offsetPartial + dimTrait];

            if (DEBUG) {
                System.err.println("\trootPrec: " + rootPrec);
            }

            if (Double.isInfinite(rootPrec)) {

                System.arraycopy(partialNodeBuffer, offsetPartial, sample, offsetSample, dimTrait);

            } else {

                final double priorPrec = partialPriorBuffer[offsetPartial + dimTrait];
                final double totalPrec = priorPrec + rootPrec;

                for (int i = 0; i < dimTrait; ++i) {
                    tmpMean[i] = (rootPrec * partialNodeBuffer[offsetPartial + i]
                            + priorPrec * partialPriorBuffer[offsetPartial + i])
                            / totalPrec;
                }

                if (DEBUG) {
                    System.err.println("\tpriorPrec: " + priorPrec);
                    System.err.println("\trootMean: " + new WrappedVector.Raw(partialNodeBuffer, offsetPartial, dimTrait));
                    System.err.println("\tprioMean: " + new WrappedVector.Raw(partialPriorBuffer, offsetPartial, dimTrait));
                    System.err.println("\tweigMean: " + new WrappedVector.Raw(tmpMean, 0, dimTrait));
                }


                final double sqrtScale = Math.sqrt(1.0 / totalPrec);

                MultivariateNormalDistribution.nextMultivariateNormalCholesky(
                        tmpMean, 0, // input mean
                        cholesky, sqrtScale, // input variance
                        sample, offsetSample, // output sample
                        tmpEpsilon);

                if (DEBUG) {
                    System.err.println("\tsample: " + new WrappedVector.Raw(sample, offsetSample, dimTrait));
                }
            }
        }

        @Override
        protected void simulateNode(final BranchNodeOperation operation, final double branchNormalization) {
            final int nodeIndex = operation.getNodeNumber();
            likelihoodDelegate.getPostOrderPartial(nodeIndex, partialNodeBuffer);

            int offsetPartial = 0;
            int offsetSample = dimNode * nodeIndex;
            int offsetParent = dimNode * operation.getParentNumber();

            final double branchPrecision = 1.0 / (operation.getBranchLength() * branchNormalization);

            if (DEBUG) {
                System.err.println("Simulate for node " + nodeIndex);
            }
             for (int trait = 0; trait < numTraits; ++trait) {

                simulateTraitForNode(nodeIndex, trait, offsetSample, offsetParent, offsetPartial, branchPrecision);

                offsetSample += dimTrait;
                offsetParent += dimTrait;
                offsetPartial += dimPartial;
            }
        }

        @Override
        protected void simulateNode(NodeOperation operation) {

        }

        protected void simulateTraitForNode(final int nodeIndex,
                                            final int traitIndex,
                                            final int offsetSample,
                                            final int offsetParent,
                                            final int offsetPartial,
                                            final double branchPrecision) {

             final double nodePrecision = partialNodeBuffer[offsetPartial + dimTrait];

             if (Double.isInfinite(nodePrecision)) {

                 if (DEBUG) {
                     System.err.println("\tCopy from node partial");
                 }

                 System.arraycopy(partialNodeBuffer, offsetPartial, sample, offsetSample, dimTrait);

             } else if (Double.isInfinite(branchPrecision)) {

                 if (DEBUG) {
                     System.err.println("\tCopy from parent sample");
                 }

                 System.arraycopy(sample, offsetParent, sample, offsetSample, dimTrait);

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

             if (DEBUG) {
                 System.err.println("\tSample value: " + new WrappedVector.Raw(sample, offsetSample, dimTrait));
             }
        }

        protected final ContinuousDataLikelihoodDelegate likelihoodDelegate;
        protected final double[] partialNodeBuffer;
        protected final double[] partialPriorBuffer;
        protected final double[] tmpMean;
    }


    class MultivariateConditionalOnTipsRealizedDelegate extends ConditionalOnTipsRealizedDelegate {

        private static final boolean DEBUG = false;

        final private PartiallyMissingInformation missingInformation;

        public MultivariateConditionalOnTipsRealizedDelegate(String name, Tree tree,
                                                             MultivariateDiffusionModel diffusionModel,
                                                             ContinuousTraitPartialsProvider dataModel,
                                                             ConjugateRootTraitPrior rootPrior,
                                                             ContinuousRateTransformation rateTransformation,
                                                             BranchRateModel rateModel,
                                                             ContinuousDataLikelihoodDelegate likelihoodDelegate) {
            super(name, tree, diffusionModel, dataModel, rootPrior, rateTransformation, rateModel, likelihoodDelegate);
            missingInformation = new PartiallyMissingInformation(tree, dataModel, likelihoodDelegate);
        }

        @Override
        protected void simulateTraitForRoot(final int offsetSample, final int offsetPartial) {

            // TODO Bad programming -- should not need to know about internal layout
            // Layout, offset, dim
            // trait, 0, dT
            // precision, dT, dT * dT
            // variance, dT + dT * dT, dT * dT
            // scalar, dT + 2 * dT * dT, 1

            // Integrate out against prior
//            final WrappedVector rootMean = new WrappedVector.Raw(partialNodeBuffer, offsetPartial, dimTrait);
//            final WrappedVector priorMean = new WrappedVector.Raw(partialPriorBuffer, offsetPartial, dimTrait);

            final DenseMatrix64F rootPrec = wrap(partialNodeBuffer, offsetPartial + dimTrait, dimTrait, dimTrait);
            final DenseMatrix64F priorPrec = new DenseMatrix64F(dimTrait, dimTrait);
            CommonOps.mult(Pd, wrap(partialPriorBuffer, offsetPartial + dimTrait, dimTrait, dimTrait), priorPrec);

            final DenseMatrix64F totalPrec = new DenseMatrix64F(dimTrait, dimTrait);
            CommonOps.add(rootPrec, priorPrec, totalPrec);

            final DenseMatrix64F totalVar = new DenseMatrix64F(dimTrait, dimTrait);
            safeInvert(totalPrec, totalVar, false);

            final double[] tmp = new double[dimTrait];
            final double[] mean = new double[dimTrait];

            for (int g = 0; g < dimTrait; ++g) {
                double sum = 0.0;
                for (int h = 0; h < dimTrait; ++h) {
                    sum += rootPrec.unsafe_get(g, h) * partialNodeBuffer[offsetPartial + h];
                    sum += priorPrec.unsafe_get(g, h) * partialPriorBuffer[offsetPartial + h];
                }
                tmp[g] = sum;
            }
            for (int g = 0; g < dimTrait; ++g) {
                double sum = 0.0;
                for (int h = 0; h < dimTrait; ++h) {
                    sum += totalVar.unsafe_get(g, h) * tmp[h];
                }
                mean[g] = sum;
            }

            final double[][] cholesky = getCholeskyOfVariance(totalVar.getData(), dimTrait);

            MultivariateNormalDistribution.nextMultivariateNormalCholesky(
                    mean, 0, // input mean
                    cholesky, 1.0, // input variance
                    sample, offsetSample, // output sample
                    tmpEpsilon);

            if (DEBUG) {
                System.err.println("Attempt to simulate root");
//                final DenseMatrix64F mean = wrap(partialNodeBuffer, offsetPartial, dimTrait, 1);
//                final DenseMatrix64F samp = wrap(sample, offsetSample, dimTrait, 1);
//                final DenseMatrix64F V = wrap(partialNodeBuffer, offsetPartial + dimTrait + dimTrait * dimTrait, dimTrait, dimTrait);

                System.err.println("Root mean: " + new WrappedVector.Raw(partialNodeBuffer, offsetPartial, dimTrait));
                System.err.println("Root prec: " + rootPrec);
                System.err.println("Priormean: " + new WrappedVector.Raw(partialPriorBuffer, offsetPartial, dimTrait));
                System.err.println("Priorprec: " + priorPrec);
                System.err.println("Totalprec: " + totalPrec);
                System.err.println("Total var: " + totalVar);


                System.err.println("draw mean: " + new WrappedVector.Raw(mean, 0, dimTrait));
//                System.err.println("V: " + totalVar);
//                System.err.println("Ch:\n" + new Matrix(cholesky));
                System.err.println("sample: " + new WrappedVector.Raw(sample, offsetSample, dimTrait));

//                System.exit(-1);
//                if (extremeValue(mean) || extremeValue(samp)) {
//                    System.exit(-1);
//                }
            }
        }

        boolean extremeValue(final DenseMatrix64F x) {
            return extremeValue(new WrappedVector.Raw(x.getData(), 0, x.getNumElements()));
        }

        boolean extremeValue(final WrappedVector x) {
            boolean valid = true;
            for (int i = 0; i < x.getDim() && valid; ++i) {
                if (Math.abs(x.get(i)) > 1E2) {
                    valid = false;
                }
            }
            return !valid;
        }

        @Override
        protected void simulateTraitForNode(final int nodeIndex,
                                            final int traitIndex,
                                            final int offsetSample,
                                             final int offsetParent,
                                             final int offsetPartial,
                                             final double branchPrecision) {

            if (nodeIndex < tree.getExternalNodeCount()) {
                simulateTraitForExternalNode(nodeIndex, traitIndex, offsetSample, offsetParent, offsetPartial, branchPrecision);
            } else {
                simulateTraitForInternalNode(offsetSample, offsetParent, offsetPartial, branchPrecision);
            }

        }

        private void simulateTraitForExternalNode(final int nodeIndex,
                                                  final int traitIndex,
                                                  final int offsetSample,
                                             final int offsetParent,
                                             final int offsetPartial,
                                             final double branchPrecision) {

            final DenseMatrix64F P0 = wrap(partialNodeBuffer, offsetPartial + dimTrait, dimTrait, dimTrait);
            final int missingCount = countFiniteDiagonals(P0);

            if (missingCount == 0) { // Completely observed

                System.arraycopy(partialNodeBuffer, offsetPartial, sample, offsetSample, dimTrait);

            } else  {

                final int zeroCount = countZeroDiagonals(P0);
                if (zeroCount == dimTrait) { //  All missing completely at random

                    final double sqrtScale = Math.sqrt(1.0 / branchPrecision);

                    MultivariateNormalDistribution.nextMultivariateNormalCholesky(
                            sample, offsetParent, // input mean
                            cholesky, sqrtScale, // input variance
                            sample, offsetSample, // output sample
                            tmpEpsilon);

                } else {

                    if (missingCount == dimTrait) { // All missing, but not completely at random

                        simulateTraitForInternalNode(offsetSample, offsetParent, offsetPartial, branchPrecision);

                    } else { // Partially observed

                        System.arraycopy(partialNodeBuffer, offsetPartial, sample, offsetSample, dimTrait); // copy observed values

                        final PartiallyMissingInformation.HashedIntArray indices = missingInformation.getMissingIndices(nodeIndex, traitIndex);
                        final int[] observed = indices.getComplement();
                        final int[] missing = indices.getArray();

                        final DenseMatrix64F V1 = new DenseMatrix64F(dimTrait, dimTrait);
                        CommonOps.scale(1.0 / branchPrecision, Vd, V1);

                        ConditionalOnPartiallyMissingTipsRealizedDelegate.ConditionalVarianceAndTranform2 transform =
                                new ConditionalOnPartiallyMissingTipsRealizedDelegate.ConditionalVarianceAndTranform2(
                                        V1, missing, observed
                                ); // TODO Cache (via delegated function)

//                        ConditionalOnPartiallyMissingTipsDelegate.ConditionalVarianceAndTranform2 transform =
//                                new ConditionalOnPartiallyMissingTipsDelegate.ConditionalVarianceAndTranform2(
//                                        Vd, missing, observed
//                                ); // TODO Cache (via delegated function)

                        final DenseMatrix64F cP0 = new DenseMatrix64F(missing.length, missing.length);
                        gatherRowsAndColumns(P0, cP0, missing, missing);

                        final WrappedVector cM2 = transform.getConditionalMean(
                                partialNodeBuffer, offsetPartial, // Tip value
                                sample, offsetParent); // Parent value

//                        final DenseMatrix64F cP1 = new DenseMatrix64F(missing.length, missing.length);
//                        CommonOps.scale(branchPrecision, transform.getConditionalPrecision(), cP1);

                        final DenseMatrix64F cP1 = transform.getConditionalPrecision();

                        final DenseMatrix64F cP2 = new DenseMatrix64F(missing.length, missing.length);
                        final DenseMatrix64F cV2 = new DenseMatrix64F(missing.length, missing.length);
                        CommonOps.add(cP0, cP1, cP2);

                        final InversionResult cc2 = safeInvert(cP2, cV2, false);
                        double[][] cC2 = getCholeskyOfVariance(cV2.getData(), missing.length);

                        MultivariateNormalDistribution.nextMultivariateNormalCholesky(
                                cM2, // input mean
                                cC2, 1.0, // input variance
                                new WrappedVector.Indexed(sample, offsetSample, missing, missing.length), // output sample
                                tmpEpsilon);


                        if (DEBUG) {
                            final WrappedVector M0 = new WrappedVector.Raw(partialNodeBuffer, offsetPartial, dimTrait);

                            final WrappedVector M1 = new WrappedVector.Raw(sample, offsetParent, dimTrait);
                            final DenseMatrix64F P1 = new DenseMatrix64F(dimTrait, dimTrait);
                            CommonOps.scale(branchPrecision, Pd, P1);

                            final WrappedVector samp = new WrappedVector.Raw(sample, offsetSample, dimTrait);

                            System.err.println("sTFEN");
                            System.err.println("M0: " + M0);
                            System.err.println("P0: " + P0);
                            System.err.println("");
                            System.err.println("M1: " + M1);
                            System.err.println("P1: " + P1);
                            System.err.println("");
//                            System.err.println("M2: " + M2);
//                            System.err.println("P2: " + P2);
//                            System.err.println("V2: " + V2);
//                            System.err.println("C2: " + new Matrix(C2));
//
//                            System.err.println("result: " + c2.getReturnCode() + " " + c2.getEffectiveDimension());
//                            System.err.println("Observed = " + new Vector(observed));
//                            System.err.println("");
                            System.err.println("");
                            System.err.println("cP0: " + cP0);
                            System.err.println("cM2: " + cM2);
                            System.err.println("cP1: " + cP1);
                            System.err.println("cP2: " + cP2);
                            System.err.println("cV2: " + cV2);
                            System.err.println("cC2: " + new Matrix(cC2));
                            System.err.println("SS: " + samp);

//                            if (extremeValue(samp)) {
//                                System.exit(-1);
//                            }

//                            System.exit(-1);

                        }
                    }
                }
            }
        }

        private void simulateTraitForInternalNode(final int offsetSample,
                                             final int offsetParent,
                                             final int offsetPartial,
                                             final double branchPrecision) {

            final WrappedVector M0 = new WrappedVector.Raw(partialNodeBuffer, offsetPartial, dimTrait);
            final DenseMatrix64F P0 = wrap(partialNodeBuffer, offsetPartial + dimTrait, dimTrait, dimTrait);

            final WrappedVector M1 = new WrappedVector.Raw(sample, offsetParent, dimTrait);
            final DenseMatrix64F P1 = new DenseMatrix64F(dimTrait, dimTrait);
            CommonOps.scale(branchPrecision, Pd, P1);

            final WrappedVector M2 = new WrappedVector.Raw(tmpMean, 0, dimTrait);
            final DenseMatrix64F P2 = new DenseMatrix64F(dimTrait, dimTrait);
            final DenseMatrix64F V2 = new DenseMatrix64F(dimTrait, dimTrait);

            CommonOps.add(P0, P1, P2);
            final InversionResult c2 = safeInvert(P2, V2, false);
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
                System.err.println("SS: " + new WrappedVector.Raw(sample, offsetSample, dimTrait));
                System.err.println("");
            }
        }
    }

    class ConditionalOnPartiallyMissingTipsRealizedDelegate extends ConditionalOnTipsRealizedDelegate {

        public ConditionalOnPartiallyMissingTipsRealizedDelegate(String name, Tree tree,
                                                         MultivariateDiffusionModel diffusionModel,
                                                         ContinuousTraitDataModel dataModel,
                                                         ConjugateRootTraitPrior rootPrior,
                                                         ContinuousRateTransformation rateTransformation,
                                                         BranchRateModel rateModel,
                                                         ContinuousDataLikelihoodDelegate likelihoodDelegate) {

            this(name, tree, diffusionModel, dataModel, rootPrior, rateTransformation, rateModel, likelihoodDelegate,
                    new PartiallyMissingInformation(tree, dataModel, likelihoodDelegate));
        }

        public ConditionalOnPartiallyMissingTipsRealizedDelegate(String name, Tree tree,
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
        protected void simulateNode(final BranchNodeOperation operation,
                                    final double branchNormalization) {
            final int nodeIndex = operation.getNodeNumber();
            likelihoodDelegate.getPostOrderPartial(nodeIndex, partialNodeBuffer);

            int offsetPartial = 0;
            int offsetSample = dimNode * nodeIndex;
            int offsetParent = dimNode * operation.getParentNumber();

            final boolean isExternal = nodeIndex < tree.getExternalNodeCount();

            final double branchPrecision = 1.0 / (operation.getBranchLength() * branchNormalization);

            for (int trait = 0; trait < numTraits; ++trait) {

                final double nodePrecision = partialNodeBuffer[offsetPartial + dimTrait];  // TODO PrecisionType.FULL

                if (!isExternal) {

                    simulateTraitForNode(nodeIndex, trait, offsetSample, offsetParent, offsetPartial, nodePrecision);

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

        public static class ConditionalVarianceAndTranform {

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
            private Matrix Sbar;
            final int[] missingIndices;
            final int[] notMissingIndices;
            final double[] tempStorage;

            final int numMissing;
            final int numNotMissing;

            private static final boolean DEBUG = false;

            public ConditionalVarianceAndTranform(final Matrix variance, final int[] missingIndices, final int[] notMissingIndices) {

                assert (missingIndices.length + notMissingIndices.length == variance.rows());
                assert (missingIndices.length + notMissingIndices.length == variance.columns());

                this.missingIndices = missingIndices;
                this.notMissingIndices = notMissingIndices;

                if (DEBUG) {
                    System.err.println("variance:\n" + variance);
                }

                Matrix S12S22Inv = null;
                Sbar = null;

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

            public double[] getConditionalMean(final double[] y, final int offsetY,
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

            public double[][] getConditionalCholesky() {
                return cholesky;
            }

            public Matrix getVariance() { return Sbar; }

            public Matrix getAffineTransform() {
                return affineTransform;
            }

            double[] getTemporageStorage() {
                return tempStorage;
            }
        }

        static class ConditionalVarianceAndTranform2 {

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

            final private DenseMatrix64F Sbar;
            final private DenseMatrix64F affineTransform;

            final int[] missingIndices;
            final int[] notMissingIndices;
            final double[] tempStorage;

            final int numMissing;
            final int numNotMissing;

            private static final boolean DEBUG = false;

            private double[][] cholesky = null;
            private DenseMatrix64F SbarInv = null;

            ConditionalVarianceAndTranform2(final DenseMatrix64F variance,
                                            final int[] missingIndices, final int[] notMissingIndices) {

                assert (missingIndices.length + notMissingIndices.length == variance.getNumRows());
                assert (missingIndices.length + notMissingIndices.length == variance.getNumCols());

                this.missingIndices = missingIndices;
                this.notMissingIndices = notMissingIndices;

                if (DEBUG) {
                    System.err.println("variance:\n" + variance);
                }

                DenseMatrix64F S22 = new DenseMatrix64F(notMissingIndices.length, notMissingIndices.length);
                gatherRowsAndColumns(variance, S22, notMissingIndices, notMissingIndices);

                if (DEBUG) {
                    System.err.println("S22:\n" + S22);
                }

                DenseMatrix64F S22Inv = new DenseMatrix64F(notMissingIndices.length, notMissingIndices.length);
                CommonOps.invert(S22, S22Inv);

                if (DEBUG) {
                    System.err.println("S22Inv:\n" + S22Inv);
                }

                DenseMatrix64F S12 = new DenseMatrix64F(missingIndices.length, notMissingIndices.length);
                gatherRowsAndColumns(variance, S12, missingIndices, notMissingIndices);

                if (DEBUG) {
                    System.err.println("S12:\n" + S12);
                }

                DenseMatrix64F S12S22Inv = new DenseMatrix64F(missingIndices.length, notMissingIndices.length);
                CommonOps.mult(S12, S22Inv, S12S22Inv);

                if (DEBUG) {
                    System.err.println("S12S22Inv:\n" + S12S22Inv);
                }

                DenseMatrix64F S12S22InvS21 = new DenseMatrix64F(missingIndices.length, missingIndices.length);
                CommonOps.multTransB(S12S22Inv, S12, S12S22InvS21);

                if (DEBUG) {
                    System.err.println("S12S22InvS21:\n" + S12S22InvS21);
                }

                Sbar = new DenseMatrix64F(missingIndices.length, missingIndices.length);
                gatherRowsAndColumns(variance, Sbar, missingIndices, missingIndices);
                CommonOps.subtract(Sbar, S12S22InvS21, Sbar);


                if (DEBUG) {
                    System.err.println("Sbar:\n" + Sbar);
                }


                this.affineTransform = S12S22Inv;
//                this.cholesky = getCholeskyOfVariance(Sbar.data, missingIndices.length);
                this.tempStorage = new double[missingIndices.length];

                this.numMissing = missingIndices.length;
                this.numNotMissing = notMissingIndices.length;

            }

            WrappedVector getConditionalMean(final double[] y, final int offsetY,
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
                        delta += affineTransform.unsafe_get(i, k) * shift[k];
                    }

                    muBar[i] = mu[offsetMu + missingIndices[i]] + delta;
                }

                return new WrappedVector.Raw(muBar, 0, numMissing);
            }

            void scatterResult(final double[] source, final int offsetSource,
                               final double[] destination, final int offsetDestination) {
                for (int i = 0; i < numMissing; ++i) {
                    destination[offsetDestination + missingIndices[i]] = source[offsetSource + i];
                }
            }

            final double[][] getConditionalCholesky() {
                if (cholesky == null) {
                    this.cholesky = getCholeskyOfVariance(Sbar.data, missingIndices.length);
                }
                return cholesky;
            }

            final DenseMatrix64F getAffineTransform() {
                return affineTransform;
            }

            final DenseMatrix64F getConditionalVariance() {
                return Sbar;
            }

            final DenseMatrix64F getConditionalPrecision() {
                if (SbarInv == null) {
                    SbarInv = new DenseMatrix64F(numMissing, numMissing);
                    CommonOps.invert(Sbar, SbarInv);
                }
                return SbarInv;
            }

            final double[] getTemporageStorage() {
                return tempStorage;
            }
        }

        public static final String PARTIAL = "partial";

        public static String getPartiallyMissingTraitName(final String traitName) {
            return PARTIAL + "." + traitName;
        }
    }

    class UnconditionalOnTipsDelegate extends AbstractRealizedContinuousTraitDelegate {

        public UnconditionalOnTipsDelegate(String name,
                                           Tree tree,
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
        protected void simulateNode(final BranchNodeOperation operation,
                                    final double branchNormalization) {
            final int nodeIndex = operation.getNodeNumber();
            int offsetSample = dimNode * nodeIndex;
            int offsetParent = dimNode * operation.getParentNumber();

            final double branchLength =  operation.getBranchLength() * branchNormalization;

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

        @Override
        protected void simulateNode(NodeOperation operation) {

        }

        private final ConjugateRootTraitPrior rootPrior;
    }
}
