/*
 * ContinuousExtensionDelegate.java
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

import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.inference.model.CompoundParameter;
import dr.math.MathUtils;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.math.matrixAlgebra.CholeskyDecomposition;
import dr.math.matrixAlgebra.WrappedVector;
import org.ejml.data.DenseMatrix64F;

/**
 * AbstractContinuousExtensionDelegate - interface for a plugin delegate for data simulation NOT on a tree.
 *
 * @author Gabriel Hassler
 */

public abstract class ContinuousExtensionDelegate {

    protected final TreeTrait treeTrait;
    protected final Tree tree;
    private final ContinuousDataLikelihoodDelegate likelihoodDelegate;
    protected final ModelExtensionProvider dataModel;

    protected final int dimTrait;
    protected final int nTaxa;

    private boolean forceResample = true;


    public ContinuousExtensionDelegate(
            ContinuousDataLikelihoodDelegate likelihoodDelegate,
            ModelExtensionProvider dataModel,
            TreeTrait treeTrait,
            Tree tree
    ) {

        this.treeTrait = treeTrait;
        this.tree = tree;
        this.likelihoodDelegate = likelihoodDelegate;
        this.dataModel = dataModel;

        this.dimTrait = dataModel.getDataDimension();
        this.nTaxa = tree.getExternalNodeCount(); //TODO: tree only used to get number of taxa
    }

    public ModelExtensionProvider getDataModel() {
        return dataModel;
    }

    public double[] getTreeTraits() {
        if (forceResample) likelihoodDelegate.fireModelChanged(); //Forces new sample
        return (double[]) treeTrait.getTrait(tree, null);
    }

    public double[] getExtendedValues() {
        double[] transformedTraits = getTransformedTraits();
        return getExtendedValues(transformedTraits);
    }

    public double[] getTransformedTraits() {
        double[] treeTraits = getTreeTraits();
        return dataModel.transformTreeTraits(treeTraits);
    }

    public double[] getExtendedValues(double[] treeValues) {
        double[] sample = new double[nTaxa * dimTrait];

        CompoundParameter dataParameter = dataModel.getParameter();
        boolean[] missingVec = dataModel.getDataMissingIndicators();

        for (int i = 0; i < nTaxa; i++) {


            IndexPartition partition = new IndexPartition(missingVec, i);
            int offset = i * dimTrait;

            for (int j : partition.obsInds) {
                int ind = j + offset;
                sample[ind] = dataParameter.getParameterValue(ind);
            }

            sampleMissingValues(sample, treeValues, partition, i);

        }

        return sample;
    }

    protected abstract void sampleMissingValues(double[] sample, double[] input, IndexPartition partition, int taxon);

    public TreeTrait getTreeTrait() {

        return treeTrait;
    }

    public Tree getTree() {
        return tree;
    }

    protected class IndexPartition {
        private final int[] obsInds;
        private final int[] misInds;
        private int nMissing;
        private int nObserved;

        private IndexPartition(boolean[] missingVector, int n) {
            int offset = n * dimTrait;

            this.nMissing = 0;

            for (int i = offset; i < offset + dimTrait; i++) {
                if (missingVector[i]) {
                    nMissing++;
                }
            }

            this.nObserved = dimTrait - nMissing;

            this.misInds = new int[nMissing];
            this.obsInds = new int[dimTrait - nMissing];
            int misInd = 0;
            int obsInd = 0;

            for (int i = offset; i < offset + dimTrait; i++) {
                if (missingVector[i]) {
                    misInds[misInd] = i - offset;
                    misInd++;
                } else {
                    obsInds[obsInd] = i - offset;
                    obsInd++;
                }
            }

        }
    }

    public static class NullExtensionDelegate extends ContinuousExtensionDelegate {

        public NullExtensionDelegate(
                ContinuousDataLikelihoodDelegate likelihoodDelegate,
                ModelExtensionProvider extensionProvider,
                TreeTrait treeTrait,
                Tree tree) {
            super(likelihoodDelegate, extensionProvider, treeTrait, tree);

        }

        @Override
        public double[] getExtendedValues() {
            return getTreeTraits();
        }

        @Override
        public double[] getExtendedValues(double[] values) {
            return values;
        }

        @Override
        protected void sampleMissingValues(double sample[], double[] input, IndexPartition partition, int taxon) {
            // do nothing
        }
    }

    public static class MultivariateNormalExtensionDelegate extends ContinuousExtensionDelegate {

        private final CompoundParameter dataParameter;

        private boolean choleskyKnown = false;
        private double[][] cholesky;
        private DenseMatrix64F extensionVariance;
        private final ModelExtensionProvider.NormalExtensionProvider dataModel;

        public MultivariateNormalExtensionDelegate(
                ContinuousDataLikelihoodDelegate likelihoodDelegate,
                TreeTrait treeTrait,
                ModelExtensionProvider.NormalExtensionProvider dataModel,
                Tree tree
        ) {
            super(likelihoodDelegate, dataModel, treeTrait, tree);

            this.dataModel = dataModel;
            this.dataParameter = dataModel.getParameter();

        }

        @Override
        public double[] getExtendedValues(double[] inputValues) {
            choleskyKnown = false;
            extensionVariance = dataModel.getExtensionVariance();
            return super.getExtendedValues(inputValues);
        }

        @Override
        protected void sampleMissingValues(double[] sample, double[] treeValues, IndexPartition partition, int taxon) {
            int offset = taxon * dimTrait;
            if (partition.nMissing == dimTrait) { // variance not diagonal, all traits missing
                double[] mean = new double[dimTrait];
                System.arraycopy(treeValues, offset, mean, 0, dimTrait);
                if (!choleskyKnown) {
                    cholesky = CholeskyDecomposition.execute(extensionVariance.getData(), 0, dimTrait);
                    choleskyKnown = true;
                }
                double[] draw = MultivariateNormalDistribution.nextMultivariateNormalCholesky(mean, cholesky);
                for (int j = offset; j < offset + dimTrait; j++) {
                    sample[j] = draw[j - offset];
                }
            } else if (partition.nMissing > 0) { // variance not diagonal, some traits missing
                ConditionalVarianceAndTransform2 transform = new ConditionalVarianceAndTransform2(extensionVariance,
                        partition.misInds, partition.obsInds);

                WrappedVector cMean = transform.getConditionalMean(
                        dataParameter.getParameter(taxon).getParameterValues(), 0, treeValues, offset);
                double[][] cChol = transform.getConditionalCholesky();

                double[] draw = MultivariateNormalDistribution.nextMultivariateNormalCholesky(cMean.getBuffer(), cChol);

                for (int j : partition.obsInds) {
                    sample[j + offset] = dataParameter.getParameterValue(j + offset);
                }
                for (int j = 0; j < partition.nMissing; j++) {
                    sample[partition.misInds[j] + offset] = draw[j];
                }
            }
        }

    }

    public static class IndependentNormalExtensionDelegate extends ContinuousExtensionDelegate {

        private final ModelExtensionProvider.NormalExtensionProvider dataModel;
        private final double[] stdev;

        public IndependentNormalExtensionDelegate(
                ContinuousDataLikelihoodDelegate likelihoodDelegate,
                TreeTrait treeTrait,
                ModelExtensionProvider.NormalExtensionProvider dataModel,
                Tree tree
        ) {
            super(likelihoodDelegate, dataModel, treeTrait, tree);

            this.dataModel = dataModel;
            this.stdev = new double[dimTrait];

        }

        @Override
        public double[] getExtendedValues(double[] inputValues) {
            DenseMatrix64F variance = dataModel.getExtensionVariance();
            for (int i = 0; i < dimTrait; i++) {
                stdev[i] = Math.sqrt(variance.get(i, i));
            }

            return super.getExtendedValues(inputValues);
        }


        @Override
        protected void sampleMissingValues(double[] sample, double[] input, IndexPartition partition, int taxon) {
            int offset = dimTrait * taxon;
            for (int j : partition.misInds) {
                int ind = j + offset;
                sample[ind] = MathUtils.nextGaussian() * stdev[j] + input[ind];
            }
        }

    }
}

