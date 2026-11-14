/*
 * BastaJNIImpl.java
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

import beagle.BeagleException;
import beagle.BeagleJNIImpl;

public class BastaJNIImpl extends BeagleJNIImpl implements BeagleBasta {

    public BastaJNIImpl(int tipCount,
                        int coalescentBufferCount,
                        int maxCoalescentIntervalCount,
                        int partialsBufferCount,
                        int compactBufferCount,
                        int stateCount,
                        int patternCount,
                        int eigenBufferCount,
                        int matrixBufferCount,
                        int categoryCount,
                        int scaleBufferCount,
                        final int[] resourceList,
                        long preferenceFlags,
                        long requirementFlags) {
        super(tipCount, partialsBufferCount, compactBufferCount, stateCount, patternCount, eigenBufferCount,
                matrixBufferCount, categoryCount, scaleBufferCount, resourceList, preferenceFlags, requirementFlags);

        int threadCount = -1;
        String tc = System.getProperty("beagle.basta.thread.count");
        if (tc != null) {
            threadCount = Integer.parseInt(tc);
        }

        allocateCoalescentBuffers(coalescentBufferCount, maxCoalescentIntervalCount, partialsBufferCount, 1, threadCount);
    }

    @Override
    public void allocateCoalescentBuffers(int coalescentBufferCount, int maxCoalescentIntervalCount, int partialsBufferCount, int initial, int threadCount) {
        int errCode = BastaJNIWrapper.INSTANCE.allocateCoalescentBuffers(instance, coalescentBufferCount,
                maxCoalescentIntervalCount, partialsBufferCount, initial, threadCount);
        if (errCode != 0) {
            throw new BeagleException("allocateCoalescentBuffers", errCode);
        }
    }

    @Override
    public void getBastaBuffer(int index, double[] buffer) {
        int errCode = BastaJNIWrapper.INSTANCE.getBastaBuffer(instance, index, buffer);
        if (errCode != 0) {
            throw new BeagleException("getCoalescentBuffers", errCode);
        }
    }

    @Override
    public void updateBastaPartialsGrad(int[] operations, int operationCount, int[] intervals, int intervalCount, int populationSizeIndex, int coalescentProbabilityIndex) {
        int errCode = BastaJNIWrapper.INSTANCE.updateBastaPartialsGrad(instance, operations, operationCount,
                intervals, intervalCount, populationSizeIndex, coalescentProbabilityIndex);
        if (errCode != 0) {
            throw new BeagleException("updateBastaPartialsGrad", errCode);
        }
    }

    @Override
    public void updateTransitionMatricesGrad(int[] transitionMatrixIndices, double[] branchLengths, int count) {
        int errCode = BastaJNIWrapper.INSTANCE.updateTransitionMatricesGrad(instance, transitionMatrixIndices, branchLengths, count);
        if (errCode != 0) {
            throw new BeagleException("updateTransitionMatricesGrad", errCode);
        }
    }

    @Override
    public void accumulateBastaPartialsGrad(int[] operations, int operationCount, int[] intervalStarts, int intervalCount, double[] intervalLengths, int populationSizeIndex, int coalescentProbabilityIndex, int eigenIndex, int partialAdjointBufferBase,
                                            int matrixAdjointBufferBase, double[] popSizeGradOut, double[] result) {
        int errCode = BastaJNIWrapper.INSTANCE.accumulateBastaPartialsGrad(instance, operations, operationCount,
                intervalStarts, intervalCount, intervalLengths, populationSizeIndex, coalescentProbabilityIndex,
                eigenIndex, partialAdjointBufferBase, matrixAdjointBufferBase, popSizeGradOut, result);
        if (errCode != 0) {
            throw new BeagleException("accumulateBastaPartialsGrad", errCode);
        }
    }

    @Override
    public void allocateCoalescentGradBuffers(int partialsCount) {
        int errCode = BastaJNIWrapper.INSTANCE.allocateCoalescentGradBuffers(instance, partialsCount);
        if (errCode != 0) {
            throw new BeagleException("allocateCoalescentGradBuffers", errCode);
        }
    }

    @Override
    public void updateBastaPartials(int[] operations, int operationCount,
                                    int[] intervals, int intervalCount,
                                    int populationSizeIndex,
                                    int coalescentProbabilityIndex) {
        int errCode = BastaJNIWrapper.INSTANCE.updateBastaPartials(instance, operations, operationCount,
                intervals, intervalCount, populationSizeIndex, coalescentProbabilityIndex);
        if (errCode != 0) {
            throw new BeagleException("updateBastaPartials", errCode);
        }
    }

    @Override
    public void accumulateBastaPartials(int[] operations, int operationCount, int[] segments, int segmentCount,
                                        double[] intervalLengths, int populationSizesIndex,
                                        int coalescentIndex, double[] result) {
        int errCode = BastaJNIWrapper.INSTANCE.accumulateBastaPartials(instance,operations, operationCount,
                segments, segmentCount, intervalLengths, populationSizesIndex, coalescentIndex, result);
        if (errCode != 0) {
            throw new BeagleException("accumulateBastaPartials", errCode);
        }
    }

    @Override
    public void accumulateEigenBasisGradient(int eigenIndex, int matrixAdjointBufferBase, double[] branchLengths,
                                              int matrixCount, boolean hasComplexEigenvalues, double[] outRateGradient) {
        int errCode = BastaJNIWrapper.INSTANCE.accumulateEigenBasisGradient(
                instance, eigenIndex, matrixAdjointBufferBase,
                branchLengths, matrixCount,
                hasComplexEigenvalues ? 1 : 0, outRateGradient);
        if (errCode != 0) {
            throw new BeagleException("accumulateEigenBasisGradient", errCode);
        }
    }


    @Override
    public void uploadBastaSlabMetadata(int[] packed, int packedLen) {
        if (packed == null) {
            throw new IllegalArgumentException("uploadBastaSlabMetadata: packed payload is null");
        }
        if (packedLen < 0 || packedLen > packed.length) {
            throw new IllegalArgumentException(
                    "uploadBastaSlabMetadata: packedLen=" + packedLen
                            + " is out of range [0, " + packed.length + "]");
        }
        int errCode = BastaJNIWrapper.INSTANCE.uploadBastaSlabMetadata(
                instance, packed, packedLen);
        if (errCode != 0) {
            throw new BeagleException("uploadBastaSlabMetadata", errCode);
        }
    }

    @Override
    public int[] getBastaSlabConstants() {
        int[] out = new int[2];
        int errCode = BastaJNIWrapper.INSTANCE.getBastaSlabConstants(instance, out);
        if (errCode == 0) {
            return out;
        }

        if (errCode == -7) {
            return null;
        }
        throw new BeagleException("getBastaSlabConstants", errCode);
    }

    public int getInstance() {
        return this.instance;
    }

}
