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
import dr.evomodel.treedatalikelihood.BufferIndexHelper;
import dr.evomodel.treedatalikelihood.continuous.cdi.ContinuousDiffusionIntegrator;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Model;

import java.util.List;

/**
 * A simple OU diffusion model delegate with branch-specific drift and constant diffusion
 * @author Marc A. Suchard
 * @author Paul Bastide
 * @version $Id$
 */
public final class OrnsteinUhlenbeckDiffusionModelDelegate extends AbstractDiffusionModelDelegate {

    // QUESTION: would it be better to nest this under DriftDiffusionModelDelegate ?
    // it was declared final, so I didn't dare to do it...

    // Here, branchRateModels represents optimal values


    private final int dim;
    private final List<BranchRateModel> branchRateModels;

    protected MatrixParameterInterface strengthOfSelectionMatrixParameter;
    private double[][] strengthOfSelectionMatrix;
    private double[][] savedStrengthOfSelectionMatrix;

//    private double[][] stationaryPrecisionMatrix;
//    private double[][] savedStationaryPrecisionMatrix;

    private final BufferIndexHelper matrixActualizationBufferHelper;

    public OrnsteinUhlenbeckDiffusionModelDelegate(Tree tree,
                                                   MultivariateDiffusionModel diffusionModel,
                                                   List<BranchRateModel> branchRateModels,
                                                   MatrixParameterInterface strengthOfSelectionMatrixParam) {
        this(tree, diffusionModel, branchRateModels, strengthOfSelectionMatrixParam, 0);
    }

    private OrnsteinUhlenbeckDiffusionModelDelegate(Tree tree,
                                                    MultivariateDiffusionModel diffusionModel,
                                                    List<BranchRateModel> branchRateModels,
                                                    MatrixParameterInterface strengthOfSelectionMatrixParam,
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
        matrixActualizationBufferHelper = new BufferIndexHelper(tree.getNodeCount(), 0, partitionNumber);
    }

    protected void calculateStrengthOfSelectionInfo(MultivariateDiffusionModel diffusionModel) {
        strengthOfSelectionMatrix = strengthOfSelectionMatrixParameter.getParameterAsMatrix();
//        stationaryPrecisionMatrix = computeStationaryVariance(strengthOfSelectionMatrix, diffusionModel.getPrecisionParameter())
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
                strengthOfSelectionMatrixParameter.getParameterValues());
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
                strengthOfSelectionMatrixParameter.getParameterValues(),
                updateCount);
    }
}