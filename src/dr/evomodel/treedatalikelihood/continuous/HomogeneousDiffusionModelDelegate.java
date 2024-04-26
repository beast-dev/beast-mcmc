/*
 * HomogenousDiffusionModelDelegate.java
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

package dr.evomodel.treedatalikelihood.continuous;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.treedatalikelihood.continuous.cdi.ContinuousDiffusionIntegrator;
import dr.math.KroneckerOperation;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

/**
 * A simple diffusion model delegate with the same diffusion model over the whole tree
 *
 * @author Marc A. Suchard
 * @author Andrew Rambaut
 * @version $Id$
 */
public final class HomogeneousDiffusionModelDelegate extends AbstractDiffusionModelDelegate {

    public HomogeneousDiffusionModelDelegate(Tree tree, MultivariateDiffusionModel diffusionModel) {
        super(tree, diffusionModel, 0);
    }

    @Override
    protected double[] getDriftRates(int[] branchIndices, int updateCount) {
        return null;
    }

    @Override
    public double[] getAccumulativeDrift(final NodeRef node, double[] priorMean, ContinuousDiffusionIntegrator cdi, int dim) {
        return priorMean;
    }

    @Override
    public double[][] getJointVariance(final double priorSampleSize, final double[][] treeVariance,
                                       final double[][] treeSharedLengths, final double[][] traitVariance) {
        return KroneckerOperation.product(treeVariance, traitVariance);
    }

    @Override
    public void getMeanTipVariances(final double priorSampleSize,
                                    final double[] treeLengths,
                                    final DenseMatrix64F traitVariance,
                                    final DenseMatrix64F varSum) {
        double sumLengths = 0;
        for (double treeLength : treeLengths) {
            sumLengths += treeLength;
        }
        sumLengths /= treeLengths.length;
        CommonOps.scale(sumLengths, traitVariance, varSum);
        CommonOps.addEquals(varSum, 1 / priorSampleSize, traitVariance);
    }
}
