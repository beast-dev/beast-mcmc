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
import dr.inference.model.Model;

import java.util.List;

/**
 * A simple diffusion model delegate with branch-specific drift and constant diffusion
 * @author Marc A. Suchard
 * @version $Id$
 */
public final class DriftDiffusionModelDelegate extends AbstractDiffusionModelDelegate {

//    private static final boolean DEBUG = true;

    private final int dim;
    private final List<BranchRateModel> branchRateModels;

    public DriftDiffusionModelDelegate(Tree tree,
                                       MultivariateDiffusionModel diffusionModel,
                                       List<BranchRateModel> branchRateModels) {
        this(tree, diffusionModel, branchRateModels, 0);
    }

    private DriftDiffusionModelDelegate(Tree tree,
                                        MultivariateDiffusionModel diffusionModel,
                                        List<BranchRateModel> branchRateModels,
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

    double[] getAccumulativeDrift(final NodeRef node) {
        final double[] drift = new double[dim];
        recursivelyAccumulateDrift(node, drift);
        return drift;
    }

    private void recursivelyAccumulateDrift(final NodeRef node, final double[] drift) {
        if (!tree.isRoot(node)) {

            final double length = tree.getBranchLength(node);

            for (int model = 0; model < dim; ++model) {
                drift[model] += branchRateModels.get(model).getBranchRate(tree, node) * length;
            }

            recursivelyAccumulateDrift(tree.getParent(node), drift);
        }
    }
}
