/*
 * DiagonalOrnsteinUhlenbeckDiffusionModelDelegate.java
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
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.inference.model.DiagonalMatrix;
import dr.inference.model.Parameter;
import org.ejml.data.DenseMatrix64F;

import java.util.List;

/**
 * A simple OU diffusion model delegate with branch-specific drift and constant diffusion
 *
 * @author Marc A. Suchard
 * @author Paul Bastide
 * @version $Id$
 */
public final class DiagonalOrnsteinUhlenbeckDiffusionModelDelegate extends AbstractOUDiffusionModelDelegate {

    // Here, branchRateModels represents optimal values

    private DiagonalMatrix strengthOfSelectionMatrixParameter;

    public DiagonalOrnsteinUhlenbeckDiffusionModelDelegate(Tree tree,
                                                           MultivariateDiffusionModel diffusionModel,
                                                           List<BranchRateModel> branchRateModels,
                                                           DiagonalMatrix strengthOfSelectionMatrixParam) {
        this(tree, diffusionModel, branchRateModels, strengthOfSelectionMatrixParam, 0);
    }

    private DiagonalOrnsteinUhlenbeckDiffusionModelDelegate(Tree tree,
                                                            MultivariateDiffusionModel diffusionModel,
                                                            List<BranchRateModel> branchRateModels,
                                                            DiagonalMatrix strengthOfSelectionMatrixParam,
                                                            int partitionNumber) {
        super(tree, diffusionModel, branchRateModels, partitionNumber);

        // Strength of selection matrix
        strengthOfSelectionMatrixParam.getDiagonalParameter().addBounds(new DiagonalMatrix.DefaultBounds(Math.log(Double.MAX_VALUE) / tree.getNodeHeight(tree.getRoot()) / 2 * Math.log(2), 0, dim));
        this.strengthOfSelectionMatrixParameter = strengthOfSelectionMatrixParam;
        addVariable(strengthOfSelectionMatrixParameter);
    }

    @Override
    public double[][] getStrengthOfSelection() {
        return strengthOfSelectionMatrixParameter.getParameterAsMatrix();
    }

    @Override
    public double[] getEigenValuesStrengthOfSelection() {
        return strengthOfSelectionMatrixParameter.getDiagonalParameter().getParameterValues();
    }

    @Override
    public double[] getEigenVectorsStrengthOfSelection() {
        return new double[0];
    }

    @Override
    public double[] getAccumulativeDrift(final NodeRef node, double[] priorMean) {
        final DenseMatrix64F drift = new DenseMatrix64F(dim, 1, true, priorMean);
        recursivelyAccumulateDrift(node, drift);
        return drift.data;
    }

    private void recursivelyAccumulateDrift(final NodeRef node, final DenseMatrix64F drift) {
        if (!tree.isRoot(node)) {

            // Compute parent
            recursivelyAccumulateDrift(tree.getParent(node), drift);


            // NOTE TO PB: Massive code duplication with work in SafeMultivariateDiagonalActualizedWithDriftIntegrator
            // Please only compute once (in SafeMultivariateDiagonalActualizedWithDriftIntegrator) and get information from cdi
            // here to accumulate

            // Actualize
            int[] branchIndice = new int[1];
            branchIndice[0] = getMatrixBufferOffsetIndex(node.getNumber());

            final double length = tree.getBranchLength(node);

            double[] actualization = new double[dim];
            computeActualizationBranch(-length, actualization);

            for (int p = 0; p < dim; ++p) {
                drift.set(p, 0, actualization[p] * drift.get(p, 0));
            }

            // Add optimal value
            double[] optVal = getDriftRates(branchIndice, 1);

            for (int p = 0; p < dim; ++p) {
                drift.set(p, 0, drift.get(p, 0) + (1 - actualization[p]) * optVal[p]);
            }

        }
    }

    private void computeActualizationBranch(double lambda, double[] C) { // NOTE TO PB: Use IntelliJ auto-formatting for consistency

        Parameter diagonals = strengthOfSelectionMatrixParameter.getDiagonalParameter(); // NOTE TO PB: avoid unnecessary copies
        for (int p = 0; p < dim; ++p) {
            C[p] = Math.exp(lambda * diagonals.getParameterValue(p));
        }
    }

    @Override
    public double[][] getJointVariance(final double priorSampleSize, final double[][] treeVariance, final double[][] treeSharedLengths, final double[][] traitVariance) {

        // Eigen of strength of selection matrix
        double[] eigVals = this.getEigenValuesStrengthOfSelection();
        int ntaxa = tree.getExternalNodeCount();
        double ti;
        double tj;
        double tij;
        double ep;
        double eq;
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
                        varTemp.set(p, q, Math.exp(-ep * ti) * Math.exp(-eq * tj) * ((Math.exp((ep + eq) * tij) - 1) / (ep + eq) + 1 / priorSampleSize) * traitVariance[p][q]);
                    }
                }
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