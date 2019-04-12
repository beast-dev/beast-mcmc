/*
 * AbstractContinuousExtensionDelegate.java
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

import dr.evolution.tree.TreeTrait;
import dr.evomodel.treedatalikelihood.continuous.ContinuousTraitPartialsProvider;
import dr.math.matrixAlgebra.WrappedVector;
import dr.math.matrixAlgebra.missingData.MissingOps;
import dr.math.matrixAlgebra.missingData.PermutationIndices;
import org.ejml.data.DenseMatrix64F;

import static dr.evomodel.treedatalikelihood.preorder.AbstractRealizedContinuousTraitDelegate.REALIZED_TIP_TRAIT;

/**
 * AbstractContinuousExtensionDelegate - interface for a plugin delegate for data simulation NOT on a tree.
 *
 * @author Gabriel Hassler
 */

public abstract class AbstractContinuousExtensionDelegate {

//    private final double[] sample;

    AbstractContinuousExtensionDelegate(ProcessSimulationDelegate.AbstractContinuousTraitDelegate
                                                treeSimulationDelegate) {


    }


    public class MultivariateNormalExtensionDelegate extends AbstractContinuousExtensionDelegate {

        private final ProcessSimulationDelegate.AbstractContinuousTraitDelegate treeSimulationDelegate;
        private final TreeTrait treeTrait;
        private final ContinuousTraitPartialsProvider partialsProvider;
        private final int dimTrait;
        private final int nTaxa;

        MultivariateNormalExtensionDelegate(
                ProcessSimulationDelegate.AbstractContinuousTraitDelegate treeSimulationDelegate,
                String traitName,
                ContinuousTraitPartialsProvider partialsProvider

        ) {
            super(treeSimulationDelegate);

            this.treeSimulationDelegate = treeSimulationDelegate;
            this.treeTrait = treeSimulationDelegate.getTreeTrait(REALIZED_TIP_TRAIT + "." + traitName);
            this.partialsProvider = partialsProvider;
            this.dimTrait = partialsProvider.getTraitDimension();
            this.nTaxa = partialsProvider.getParameter().getParameterCount();

        }

        public double[] getExtendedValues() {


            double[] treeValues = (double[]) treeTrait.getTrait(treeSimulationDelegate.tree, null); //TODO: make sure modifying this array doesn't break anything

            for (int i = 0; i < nTaxa; i++) {

                double[] tipPartial = partialsProvider.getTipPartial(i, false); //TODO: should fullyObserved always be false?
                DenseMatrix64F V = MissingOps.wrap(tipPartial, dimTrait + dimTrait * dimTrait, dimTrait, dimTrait);
                final PermutationIndices inds = new PermutationIndices(V);
                int[] missingInds = inds.getInfiniteIndices();
                int[] obsInds = inds.getNonZeroFiniteIndices();
                if (inds.getNumberOfInfiniteDiagonals() + inds.getNumberOfNonZeroFiniteDiagonals() != dimTrait) {
                    throw new RuntimeException("Not currently implemented for cases where traits are both latent and fully observed.");
                    //TODO: construct obsInds such that it can handle this issue use case.
                }

                ConditionalVarianceAndTransform2 transform = new ConditionalVarianceAndTransform2(V, missingInds, obsInds);

                int offset = i * dimTrait;

                WrappedVector cMean = transform.getConditionalMean(tipPartial, 0, treeValues, offset);
                double[][] cChol = transform.getConditionalCholesky();


            }


            return null;
        }
    }

}

