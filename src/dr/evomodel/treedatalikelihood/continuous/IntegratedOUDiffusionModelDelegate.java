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

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.continuous.MultivariateElasticModel;
import dr.evomodel.treedatalikelihood.continuous.cdi.ContinuousDiffusionIntegrator;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import java.util.List;

import static dr.math.matrixAlgebra.missingData.MissingOps.wrap;

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
    public DenseMatrix64F getGradientVarianceWrtVariance(NodeRef node,
                                                         ContinuousDiffusionIntegrator cdi,
                                                         ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                                         DenseMatrix64F gradient) {
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
    public double[][] getJointVariance(final double priorSampleSize,
                                       final double[][] treeVariance, final double[][] treeSharedLengths,
                                       final double[][] traitVariance) {

        double[] eigVals = this.getEigenValuesStrengthOfSelection();
        DenseMatrix64F V = wrap(this.getEigenVectorsStrengthOfSelection(), 0, dim, dim);
        DenseMatrix64F Vinv = new DenseMatrix64F(dim, dim);
        CommonOps.invert(V, Vinv);

        DenseMatrix64F transTraitVariance = new DenseMatrix64F(traitVariance);

        DenseMatrix64F tmp = new DenseMatrix64F(dim, dim);
        CommonOps.mult(Vinv, transTraitVariance, tmp);
        CommonOps.multTransB(tmp, Vinv, transTraitVariance);

        // Computation of matrix
        int ntaxa = tree.getExternalNodeCount();
        double ti;
        double tj;
        double tij;
        double ep;
        double eq;
        double var;
        DenseMatrix64F varTemp = new DenseMatrix64F(dim, dim);
        double[][] jointVariance = new double[dim * ntaxa][dim * ntaxa];
        for (int i = 0; i < ntaxa; ++i) {
            for (int j = 0; j < ntaxa; ++j) {
                ti = treeSharedLengths[i][i];
                tj = treeSharedLengths[j][j];
                tij = treeSharedLengths[i][j];
                for (int p = 0; p < dim; ++p) {
                    for (int q = 0; q < dim; ++q) {
                        ep = eigVals[p];
                        eq = eigVals[q];
                        var = tij / ep / eq;
                        var += (1 - Math.exp(ep * tij)) * Math.exp(-ep * ti) / ep / ep / eq;
                        var += (1 - Math.exp(eq * tij)) * Math.exp(-eq * tj) / ep / eq / eq;
                        var -= (1 - Math.exp((ep + eq) * tij)) * Math.exp(-ep * ti) * Math.exp(-eq * tj) / ep / eq / (ep + eq);
                        var += (1 - Math.exp(-ep * ti)) * (1 - Math.exp(-eq * tj)) / ep / eq / priorSampleSize;
                        var += 1 / priorSampleSize;
                        varTemp.set(p, q, var * transTraitVariance.get(p, q));
                    }
                }
                CommonOps.mult(V, varTemp, tmp);
                CommonOps.multTransB(tmp, V, varTemp);
                for (int p = 0; p < dim; ++p) {
                    for (int q = 0; q < dim; ++q) {
                        jointVariance[i * dim + p][j * dim + q] = varTemp.get(p, q);
                    }
                }
            }
        }
        return jointVariance;
    }
}
