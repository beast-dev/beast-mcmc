/*
 * DIffusionProcessDelegate.java
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
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.treedatalikelihood.continuous.cdi.ContinuousDiffusionIntegrator;
import dr.inference.model.Model;
import org.ejml.data.DenseMatrix64F;

/**
 * Implementations of this interface are used to delegate control of diffusion
 * models along branches of a tree.
 *
 * @author Marc Suchard
 * @author Andrew Rambaut
 * @version $Id$
 */
public interface DiffusionProcessDelegate extends Model {

    int getDiffusionModelCount();

    int getEigenBufferCount();

    int getEigenBufferOffsetIndex(int i);

    int getMatrixBufferCount();

    public int getMatrixBufferOffsetIndex(int i);

    void flipMatrixBufferOffset(int i);

    MultivariateDiffusionModel getDiffusionModel(int index);

    int getMatrixIndex(int branchIndex);

    void setDiffusionModels(ContinuousDiffusionIntegrator cdi, boolean flipBuffers);

    void updateDiffusionMatrices(ContinuousDiffusionIntegrator cdi, int[] branchIndices, double[] edgeLengths,
                                 int updateCount, boolean flipBuffers);

    boolean hasDrift();

    boolean hasActualization();

    boolean hasDiagonalActualization();

    boolean isIntegratedProcess();

    DenseMatrix64F getGradientVarianceWrtVariance(NodeRef node, ContinuousDiffusionIntegrator cdi, ContinuousDataLikelihoodDelegate likelihoodDelegate, DenseMatrix64F gradient);

    double[] getGradientDisplacementWrtRoot(NodeRef node, ContinuousDiffusionIntegrator cdi, ContinuousDataLikelihoodDelegate likelihoodDelegate, DenseMatrix64F gradient);

    void storeState();

    void restoreState();

    double[] getAccumulativeDrift(final NodeRef node, double[] priorMean, ContinuousDiffusionIntegrator cdi, int dim);

    double[][] getJointVariance(final double priorSampleSize, final double[][] treeVariance, final double[][] treeSharedLengths, final double[][] traitVariance);

    void getMeanTipVariances(final double priorSampleSize, final double[] treeLengths, final DenseMatrix64F traitVariance, final DenseMatrix64F varSum);
}
