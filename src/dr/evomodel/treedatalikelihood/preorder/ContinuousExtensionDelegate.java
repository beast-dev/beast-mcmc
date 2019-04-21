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
import dr.math.distributions.MultivariateNormalDistribution;
import dr.math.matrixAlgebra.WrappedVector;
import org.ejml.data.DenseMatrix64F;

/**
 * AbstractContinuousExtensionDelegate - interface for a plugin delegate for data simulation NOT on a tree.
 *
 * @author Gabriel Hassler
 */

public class ContinuousExtensionDelegate {

    protected final TreeTrait treeTrait;
    protected final Tree tree;
    private final ContinuousDataLikelihoodDelegate likelihoodDelegate;

    public ContinuousExtensionDelegate(
            ContinuousDataLikelihoodDelegate likelihoodDelegate,
            TreeTrait treeTrait,
            Tree tree
    ) {

        this.treeTrait = treeTrait;
        this.tree = tree;
        this.likelihoodDelegate = likelihoodDelegate;
    }

    public double[] getExtendedValues() {
        likelihoodDelegate.fireModelChanged(); //Forces new sample

        return (double[]) treeTrait.getTrait(tree, null);
    }

    public double[] getExtendedValues(double[] tipTraits) {
        return tipTraits;
    }

    public TreeTrait getTreeTrait() {

        return treeTrait;
    }

    public Tree getTree() {
        return tree;
    }


    public static class MultivariateNormalExtensionDelegate extends ContinuousExtensionDelegate {
        private final double[] sample;
        private final ModelExtensionProvider.NormalExtensionProvider dataModel;
        private final int dimTrait;
        private final int nTaxa;

        public MultivariateNormalExtensionDelegate(
                ContinuousDataLikelihoodDelegate likelihoodDelegate,
                TreeTrait treeTrait,
                ModelExtensionProvider.NormalExtensionProvider dataModel,
                Tree tree
        ) {
            super(likelihoodDelegate, treeTrait, tree);
            this.dataModel = dataModel;
            this.dimTrait = dataModel.getTraitDimension();
            this.nTaxa = tree.getExternalNodeCount();
            this.sample = new double[nTaxa * dimTrait];

        }

        @Override
        public double[] getExtendedValues() {
            double[] treeValues = super.getExtendedValues();
            return getExtendedValues(treeValues);
        }

        @Override
        public double[] getExtendedValues(double[] treeValues) {

            CompoundParameter dataParameter = dataModel.getParameter();
            DenseMatrix64F extensionVar = dataModel.getExtensionVariance();
            boolean[] missingVec = dataModel.getMissingVector();

            int offset = 0;


            for (int i = 0; i < nTaxa; i++) {

                IndexPartition partition = new IndexPartition(missingVec, i);

                if (partition.nObserved == dimTrait) {
                    for (int j = offset; j < offset + dimTrait; j++) {
                        sample[j] = dataParameter.getParameterValue(j);
                    }
                } else if (partition.nMissing == dimTrait) {
                    //TODO: completely missing
                    throw new RuntimeException("Not currently implemented for all traits missing.");
                } else {
                    ConditionalVarianceAndTransform2 transform = new ConditionalVarianceAndTransform2(extensionVar,
                            partition.misInds, partition.obsInds);

                    WrappedVector cMean = transform.getConditionalMean(
                            dataParameter.getParameter(i).getParameterValues(), 0, treeValues, offset);
                    double[][] cChol = transform.getConditionalCholesky();

                    double[] draw = MultivariateNormalDistribution.nextMultivariateNormalCholesky(cMean.getBuffer(), cChol);

                    for (int j : partition.obsInds) {
                        sample[j + offset] = dataParameter.getParameterValue(j + offset);
                    }
                    for (int j = 0; j < partition.nMissing; j++) {
                        sample[partition.misInds[j] + offset] = draw[j];
                    }
                }

                offset += dimTrait;


            }


            return sample;
        }

        private class IndexPartition {
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

    }

}

