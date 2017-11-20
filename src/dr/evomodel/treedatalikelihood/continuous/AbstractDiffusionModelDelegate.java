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

import dr.evolution.tree.Tree;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.treedatalikelihood.BufferIndexHelper;
import dr.evomodel.treedatalikelihood.continuous.cdi.ContinuousDiffusionIntegrator;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

import java.io.Serializable;

/**
 * A simple diffusion model delegate with the same diffusion model over the whole tree
 * @author Marc A. Suchard
 * @author Andrew Rambaut
 * @version $Id$
 */
public abstract class AbstractDiffusionModelDelegate extends AbstractModel implements DiffusionProcessDelegate, Serializable {

//    private static final boolean DEBUG = false;

    final Tree tree;
    private final MultivariateDiffusionModel diffusionModel;

    private final BufferIndexHelper eigenBufferHelper;
    private final BufferIndexHelper matrixBufferHelper;

    AbstractDiffusionModelDelegate(Tree tree, MultivariateDiffusionModel diffusionModel,
                                          int partitionNumber) {

        super("AbstractDiffusionModelDelegate");

        this.tree = tree;
        this.diffusionModel = diffusionModel;
        addModel(diffusionModel);

        // two eigen buffers for each decomposition for store and restore.
        eigenBufferHelper = new BufferIndexHelper(1, 0, partitionNumber);

        // two matrices for each node less the root
        matrixBufferHelper = new BufferIndexHelper(tree.getNodeCount(), 0, partitionNumber);
    }

    @Override
    public int getEigenBufferCount() {
        return eigenBufferHelper.getBufferCount();
    }

    @Override
    public int getMatrixBufferCount() {
        return matrixBufferHelper.getBufferCount();
    }

    @Override
    public int getDiffusionModelCount() {
        return 1;
    }

    @Override
    public MultivariateDiffusionModel getDiffusionModel(int index) {
        assert(index == 0);
        return diffusionModel;
    }

    @Override
    public int getMatrixIndex(int branchIndex) {
        return matrixBufferHelper.getOffsetIndex(branchIndex);
    }

    @Override
    public void setDiffusionModels(ContinuousDiffusionIntegrator cdi, boolean flip) {
        if (flip) {
            eigenBufferHelper.flipOffset(0);
        }

        cdi.setDiffusionPrecision(eigenBufferHelper.getOffsetIndex(0),
                diffusionModel.getPrecisionParameter().getParameterValues(),
                Math.log(diffusionModel.getDeterminantPrecisionMatrix())

        );
    }

    @Override
    public void updateDiffusionMatrices(ContinuousDiffusionIntegrator cdi, int[] branchIndices, double[] edgeLengths,
                                        int updateCount, boolean flip) {

        int[] probabilityIndices = new int[updateCount];

        for (int i = 0; i < updateCount; i++) {
            if (flip) {
                matrixBufferHelper.flipOffset(branchIndices[i]);
            }
            probabilityIndices[i] = matrixBufferHelper.getOffsetIndex(branchIndices[i]);
        }

        cdi.updateBrownianDiffusionMatrices(
                eigenBufferHelper.getOffsetIndex(0),
                probabilityIndices,
                edgeLengths,
                getDriftRates(branchIndices, updateCount),
                updateCount);
    }

    protected abstract double[] getDriftRates(int[] branchIndices, int updateCount);

    @Override
    public boolean hasDrift() { return false; }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == diffusionModel) {
            fireModelChanged(model);
        } else {
            throw new RuntimeException("Unknown model");
        }
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {

    }

    @Override
    public void storeState() {
        eigenBufferHelper.storeState();
        matrixBufferHelper.storeState();
    }

    @Override
    public void restoreState() {
        eigenBufferHelper.restoreState();
        matrixBufferHelper.restoreState();
    }

    @Override
    protected void acceptState() {

    }
}
