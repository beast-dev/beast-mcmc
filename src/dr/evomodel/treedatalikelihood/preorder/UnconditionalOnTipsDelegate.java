/*
 * UnconditionalOnTipsDelegate.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.evomodel.treedatalikelihood.preorder;

import dr.evolution.tree.MutableTreeModel;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.treedatalikelihood.continuous.ConjugateRootTraitPrior;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.continuous.ContinuousRateTransformation;
import dr.evomodel.treedatalikelihood.continuous.ContinuousTraitDataModel;
import dr.math.distributions.MultivariateNormalDistribution;

/**
 * @author Marc A. Suchard
 */
public class UnconditionalOnTipsDelegate extends AbstractRealizedContinuousTraitDelegate {

    public UnconditionalOnTipsDelegate(String name,
                                       MutableTreeModel tree,
                                       MultivariateDiffusionModel diffusionModel,
                                       ContinuousTraitDataModel dataModel,
                                       ConjugateRootTraitPrior rootPrior,
                                       ContinuousRateTransformation rateTransformation,
                                       ContinuousDataLikelihoodDelegate likelihoodDelegate) {
        super(name, tree, diffusionModel, dataModel, rootPrior, rateTransformation, likelihoodDelegate);

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

    protected void simulateNode(final BranchNodeOperation operation,
                                final double branchNormalization) {
        final int nodeIndex = operation.getNodeNumber();
        int offsetSample = dimNode * nodeIndex;
        int offsetParent = dimNode * operation.getParentNumber();

        final double branchLength = operation.getBranchLength() * branchNormalization;

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
    protected void simulateNode(final int parentNumber,
                                         final int nodeNumber,
                                         final int nodeMatrix,
                                         final int siblingNumber,
                                         final int siblingMatrix) {
        throw new RuntimeException("Not yet implemented -- see above");
    }

    private final ConjugateRootTraitPrior rootPrior;
}
