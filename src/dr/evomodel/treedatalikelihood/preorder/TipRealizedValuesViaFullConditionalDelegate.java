/*
 * TipRealizedValuesViaFullConditionalDelegate.java
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

import dr.evolution.tree.Tree;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.treedatalikelihood.continuous.ConjugateRootTraitPrior;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.continuous.ContinuousRateTransformation;
import dr.evomodel.treedatalikelihood.continuous.ContinuousTraitPartialsProvider;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.math.matrixAlgebra.WrappedMatrix;
import dr.math.matrixAlgebra.WrappedVector;

/**
 * @author Marc A. Suchard
 */
public class TipRealizedValuesViaFullConditionalDelegate extends AbstractValuesViaFullConditionalDelegate {

    public TipRealizedValuesViaFullConditionalDelegate(String name, Tree tree,
                                                       MultivariateDiffusionModel diffusionModel,
                                                       ContinuousTraitPartialsProvider dataModel,
                                                       ConjugateRootTraitPrior rootPrior,
                                                       ContinuousRateTransformation rateTransformation,
//                                                       BranchRateModel rateModel,
                                                       ContinuousDataLikelihoodDelegate likelihoodDelegate) {
        super(name, tree, diffusionModel, dataModel, rootPrior, rateTransformation, likelihoodDelegate);
    }

    public static String getName(String name) {
        return "tipSample." + name;
    }

    public String getTraitName(String name) {
        return getName(name);
    }

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
                new WrappedMatrix.ArrayOfArray(cholesky), 1.0, // input variance
                output, // output sample
                buffer);
    }
}
