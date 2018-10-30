/*
 * IntegratedOUDiffusionModelDelegate.java
 *
 * Copyright (c) 2002-2018 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.continuous.MultivariateElasticModel;
import dr.evomodel.treedatalikelihood.continuous.cdi.ContinuousDiffusionIntegrator;
import org.ejml.data.DenseMatrix64F;

import java.util.List;

/**
 * A simple OU diffusion model delegate with branch-specific drift and constant diffusion
 *
 * @author Marc A. Suchard
 * @author Paul Bastide
 * @version $Id$
 */
public class IntegratedOUDiffusionModelDelegate extends OUDiffusionModelDelegate {

    public IntegratedOUDiffusionModelDelegate(Tree tree,
                                              MultivariateDiffusionModel diffusionModel,
                                              List<BranchRateModel> branchRateModels,
                                              MultivariateElasticModel elasticModel) {
        super(tree, diffusionModel, branchRateModels, elasticModel);
    }

    @Override
    public boolean hasDrift() {
        return true;
    }

    @Override
    public boolean hasActualization() {
        return true;
    }

    @Override
    public boolean hasDiagonalActualization() {
        return false;
    }

    @Override
    public boolean isIntegratedProcess() {
        return true;
    }

    public boolean isSymmetric() {
        return false;
    }

    @Override
    public void updateDiffusionMatrices(ContinuousDiffusionIntegrator cdi, int[] branchIndices, double[] edgeLengths,
                                        int updateCount, boolean flip) {

        int[] probabilityIndices = new int[updateCount];

        for (int i = 0; i < updateCount; i++) {
            if (flip) {
                flipMatrixBufferOffset(branchIndices[i]);
            }
            probabilityIndices[i] = getMatrixBufferOffsetIndex(branchIndices[i]);
        }

        cdi.updateIntegratedOrnsteinUhlenbeckDiffusionMatrices(
                getEigenBufferOffsetIndex(0),
                probabilityIndices,
                edgeLengths,
                getDriftRates(branchIndices, updateCount),
                getEigenValuesStrengthOfSelection(),
                getEigenVectorsStrengthOfSelection(),
                updateCount);
    }

    @Override
    public void getGradientPrecision(double scalar, DenseMatrix64F gradient) {
        throw new RuntimeException("not yet implemented");
    }

    @Override
    public double[][] getJointVariance(final double priorSampleSize,
                                       final double[][] treeVariance, final double[][] treeSharedLengths,
                                       final double[][] traitVariance) {
        throw new RuntimeException("Not yet implemented");
    }

    private double[][] getJointVarianceFull(final double priorSampleSize,
                                            final double[][] treeVariance, final double[][] treeSharedLengths,
                                            final double[][] traitVariance) {
        throw new RuntimeException("Not yet implemented");
    }

    private double[][] getJointVarianceDiagonal(final double priorSampleSize,
                                                final double[][] treeVariance, final double[][] treeSharedLengths,
                                                final double[][] traitVariance) {
        throw new RuntimeException("Not yet implemented");
    }
}
