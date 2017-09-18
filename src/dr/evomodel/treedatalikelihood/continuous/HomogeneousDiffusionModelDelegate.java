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
import dr.evomodel.treedatalikelihood.continuous.cdi.ContinuousDiffusionIntegrator;

/**
 * A simple diffusion model delegate with the same diffusion model over the whole tree
 * @author Marc A. Suchard
 * @author Andrew Rambaut
 * @version $Id$
 */
public final class HomogeneousDiffusionModelDelegate extends AbstractDiffusionModelDelegate {

//    private static final boolean DEBUG = false;

    public HomogeneousDiffusionModelDelegate(Tree tree, MultivariateDiffusionModel diffusionModel) {
        this(tree, diffusionModel, 0);
    }

    private HomogeneousDiffusionModelDelegate(Tree tree, MultivariateDiffusionModel diffusionModel,
                                              int partitionNumber) {
        super(tree, diffusionModel, partitionNumber);
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
                edgeLengths, null,
                updateCount);
    }
}
