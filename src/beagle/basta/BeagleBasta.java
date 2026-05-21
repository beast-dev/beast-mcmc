/*
 * BeagleBasta.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package beagle.basta;

import beagle.Beagle;

public interface BeagleBasta extends Beagle {

    int BASTA_OPERATION_SIZE = 8;

    void updateBastaPartials(final int[] operations,
                             int operationCount,
                             final int[] intervals,
                             int intervalCount,
                             int populationSizeIndex,
                             int coalescentProbabilityIndex);

    void accumulateBastaPartials(final int[] operations,
                                 int operationCount,
                                 final int[] intervalStarts,
                                 int intervalCount,
                                 final double[] intervalLengths,
                                 int populationSizesIndex,
                                 int coalescentIndex,
                                 double[] result);

    void allocateCoalescentBuffers(int coalescentBufferCount, int maxCoalescentIntervalCount, int partialsBufferCount, int initial, int threadCount);

    void getBastaBuffer(int index, double[] buffer);

    void updateBastaPartialsGrad(int[] operations, int operationCount, int[] intervals, int intervalCount, int populationSizeIndex, int coalescentProbabilityIndex);

    void updateTransitionMatricesGrad(int[] transitionMatrixIndices, double[] branchLengths, int count);

    void accumulateBastaPartialsGrad(int[] operations, int operationCount,  final int[] intervalStarts, int intervalCount, final double[] intervalLengths, int populationSizeIndex, int coalescentProbabilityIndex, double[] result);

    void allocateCoalescentGradBuffers(int partialsCount);

    void updateBastaPartialsPopSizeGrad(int[] operations, int operationCount, int[] intervals, int intervalCount, int populationSizeIndex, int coalescentProbabilityIndex);

    void accumulateBastaPartialsPopSizeGrad(int[] operations, int operationCount, final int[] intervalStarts, int intervalCount, final double[] intervalLengths, int populationSizeIndex, int coalescentProbabilityIndex, double[] result);

    void getMatrixAdjoint(int matrixIndex, double[] buffer);

    void getPopulationSizeGradient(double[] buffer);

    void setExpmKernels(double[] kernels);

    void accumulateExpmGradient(double[] out);

    void transformMatrixAdjoints(int matrixCount, double[] out);

    void backTransformEigenBasisGradient(double[] eigenBasisGrad, double[] out);

    void accumulateEigenBasisGradient(double[] eigenValues, double[] branchLengths,
                                      int matrixCount, boolean hasComplexEigenvalues,
                                      double[] outRateGradient);


    void uploadBastaSlabMetadata(int[] packed, int packedLen);

    int[] getBastaSlabConstants();
}
