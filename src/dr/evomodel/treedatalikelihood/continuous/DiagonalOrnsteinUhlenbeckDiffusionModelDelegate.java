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
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.treedatalikelihood.continuous.cdi.ContinuousDiffusionIntegrator;
import dr.inference.model.DiagonalMatrix;
import dr.inference.model.Model;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import java.util.List;

/**
 * A simple OU diffusion model delegate with branch-specific drift and constant diffusion
 * @author Marc A. Suchard
 * @author Paul Bastide
 * @version $Id$
 */
public final class DiagonalOrnsteinUhlenbeckDiffusionModelDelegate extends AbstractDiffusionModelDelegate {

    // Here, branchRateModels represents optimal values


    private final int dim;
    private final List<BranchRateModel> branchRateModels;

    //    protected MatrixParameterInterface strengthOfSelectionMatrixParameter;
    protected DiagonalMatrix strengthOfSelectionMatrixParameter;
    private double[][] strengthOfSelectionMatrix;
//    private double[][] savedStrengthOfSelectionMatrix;

//    private double[][] stationaryPrecisionMatrix;
//    private double[][] savedStationaryPrecisionMatrix;

//    private final BufferIndexHelper matrixActualizationBufferHelper;

    public DiagonalOrnsteinUhlenbeckDiffusionModelDelegate(Tree tree,
                                                   MultivariateDiffusionModel diffusionModel,
                                                   List<BranchRateModel> branchRateModels,
//                                                   MatrixParameterInterface strengthOfSelectionMatrixParam,
                                                   DiagonalMatrix strengthOfSelectionMatrixParam) {
        this(tree, diffusionModel, branchRateModels, strengthOfSelectionMatrixParam, 0);
    }

    private DiagonalOrnsteinUhlenbeckDiffusionModelDelegate(Tree tree,
                                                    MultivariateDiffusionModel diffusionModel,
                                                    List<BranchRateModel> branchRateModels,
//                                                    MatrixParameterInterface strengthOfSelectionMatrixParam,
                                                    DiagonalMatrix strengthOfSelectionMatrixParam,
                                                    int partitionNumber) {
        super(tree, diffusionModel, partitionNumber);
        this.branchRateModels = branchRateModels;

        dim = diffusionModel.getPrecisionParameter().getColumnDimension();

        if (branchRateModels != null) {

            for (BranchRateModel rateModel : branchRateModels) {
                addModel(rateModel);
            }

            if (branchRateModels.size() != dim) {
                throw new IllegalArgumentException("Invalid dimensions");
            }
        }

        // Strength of selection matrix
        this.strengthOfSelectionMatrixParameter = strengthOfSelectionMatrixParam;
        calculateStrengthOfSelectionInfo(diffusionModel);
        addVariable(strengthOfSelectionMatrixParameter);

        // two more matrices for each node less the root
//        matrixActualizationBufferHelper = new BufferIndexHelper(tree.getNodeCount(), 0, partitionNumber);
    }

    protected void calculateStrengthOfSelectionInfo(MultivariateDiffusionModel diffusionModel) {
        strengthOfSelectionMatrix = strengthOfSelectionMatrixParameter.getParameterAsMatrix();
//        stationaryPrecisionMatrix = computeStationaryVariance(strengthOfSelectionMatrix, diffusionModel.getPrecisionParameter())
    }

    public double[] getStrengthOfSelection() {
        return strengthOfSelectionMatrixParameter.getDiagonalParameter().getParameterValues();
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {

        if (branchRateModels.contains(model)) {
            fireModelChanged(model);
        } else {
            super.handleModelChangedEvent(model, object, index);
        }
    }

    @Override
    public boolean hasDrift() { return true; }

    @Override
    public boolean hasActualization() { return true; }

    @Override
    protected double[] getDriftRates(int[] branchIndices, int updateCount) {

        final double[] drift = new double[updateCount * dim];  // TODO Reuse?

        if (branchRateModels != null) {

            int offset = 0;
            for (int i = 0; i < updateCount; ++i) {

                final NodeRef node = tree.getNode(branchIndices[i]); // TODO Check if correct node

                for (int model = 0; model < dim; ++model) {
                    drift[offset] = branchRateModels.get(model).getBranchRate(tree, node);
                    ++offset;
                }
            }
        }
        return drift;
    }

    @Override
    public void setDiffusionModels(ContinuousDiffusionIntegrator cdi, boolean flip) {
        super.setDiffusionModels(cdi, flip);

        cdi.setDiffusionStationaryVariance(getEigenBufferOffsetIndex(0),
                strengthOfSelectionMatrixParameter.getDiagonalParameter().getParameterValues());
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

        cdi.updateOrnsteinUhlenbeckDiffusionMatrices(
                getEigenBufferOffsetIndex(0),
                probabilityIndices,
                edgeLengths,
                getDriftRates(branchIndices, updateCount),
                strengthOfSelectionMatrixParameter.getDiagonalParameter().getParameterValues(),
                updateCount);
    }

    double[] getAccumulativeDrift(final NodeRef node, double[] priorMean) {
        final DenseMatrix64F drift = new DenseMatrix64F(dim, 1, true, priorMean);
        recursivelyAccumulateDrift(node, drift);
        return drift.data;
    }

    private void recursivelyAccumulateDrift(final NodeRef node, final DenseMatrix64F drift) {
        if (!tree.isRoot(node)) {

            // Compute parent
            recursivelyAccumulateDrift(tree.getParent(node), drift);

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
                drift.set(p, 0, drift.get(p, 0) + (1 - actualization[p]) *  optVal[p]);
            }

        }
    }

    private void computeActualizationBranch(double lambda, double[] C){
        double[] A = strengthOfSelectionMatrixParameter.getDiagonalParameter().getParameterValues();
        for (int p = 0; p < dim; ++p) {
            C[p] = Math.exp(lambda * A[p]);
        }
    }
}