/*
 * IntegratedBMDiffusionModelDelegate.java
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

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.treedatalikelihood.continuous.cdi.ContinuousDiffusionIntegrator;
import dr.math.KroneckerOperation;
import org.ejml.data.DenseMatrix64F;

import java.util.List;

/**
 * A simple OU diffusion model delegate with branch-specific drift and constant diffusion
 *
 * @author Marc A. Suchard
 * @author Paul Bastide
 * @version $Id$
 */
public class IntegratedBMDiffusionModelDelegate extends AbstractDriftDiffusionModelDelegate {

    public IntegratedBMDiffusionModelDelegate(Tree tree,
                                              MultivariateDiffusionModel diffusionModel,
                                              List<BranchRateModel> branchRateModels) {
        this(tree, diffusionModel, branchRateModels, 0);
    }

    private IntegratedBMDiffusionModelDelegate(Tree tree,
                                               MultivariateDiffusionModel diffusionModel,
                                               List<BranchRateModel> branchRateModels,
                                               int partitionNumber) {
        super(tree, diffusionModel, branchRateModels, partitionNumber);
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

        cdi.updateIntegratedBrownianMotionDiffusionMatrices(
                getEigenBufferOffsetIndex(0),
                probabilityIndices,
                edgeLengths,
                getDriftRates(branchIndices, updateCount),
                updateCount);
    }

    @Override
    public void getGradientPrecision(double scalar, DenseMatrix64F gradient) {
        throw new RuntimeException("not yet implemented");
    }

    @Override
    public double[] getAccumulativeDrift(final NodeRef node, double[] priorMean, ContinuousDiffusionIntegrator cdi, int dim) {
        double[] driftFull = super.getAccumulativeDrift(node, priorMean, cdi, dim);
        assert dim % 2 == 0 : "dimTrait should be twice dimProcess.";
        int dimProcess = dim / 2;
        double[] drift = new double[dimProcess];
        System.arraycopy(driftFull, dimProcess, drift, 0, dimProcess);
        return drift;
    }

    @Override
    public double[][] getJointVariance(final double priorSampleSize, final double[][] treeVariance,
                                       final double[][] treeSharedLengths, final double[][] traitVariance) {
        double[][] integratedTreeVariance = computeIntegratedTreeVariance(treeSharedLengths, priorSampleSize);
        return KroneckerOperation.product(integratedTreeVariance, traitVariance);
    }

    private double[][] computeIntegratedTreeVariance(double[][] treeSharedLengths, double priorSampleSize) {
        int ntaxa = treeSharedLengths.length;
        double[][] integratedTreeVariance = new double[ntaxa][ntaxa];
        for (int i = 0; i < ntaxa; i++) {
            double ti = treeSharedLengths[i][i];
            for (int j = 0; j < ntaxa; j++){
                double tj = treeSharedLengths[j][j];
                double tij = treeSharedLengths[i][j];
                integratedTreeVariance[i][j] = (1 + ti * tj) / priorSampleSize + tij * (ti * tj + tij * (tij / 3 - (ti + tj) / 2));
            }
        }
        return integratedTreeVariance;
    }
}
