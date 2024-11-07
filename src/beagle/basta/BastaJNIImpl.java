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

        allocateCoalescentBuffers(coalescentBufferCount, maxCoalescentIntervalCount, partialsBufferCount,1);
    }

    @Override
    public void allocateCoalescentBuffers(int coalescentBufferCount, int maxCoalescentIntervalCount, int partialsBufferCount, int initial) {
        int errCode = BastaJNIWrapper.INSTANCE.allocateCoalescentBuffers(instance, coalescentBufferCount,
                maxCoalescentIntervalCount, partialsBufferCount, initial);
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
    public void accumulateBastaPartialsGrad(int[] operations, int operationCount, int[] intervalStarts, int intervalCount, double[] intervalLengths, int populationSizeIndex, int coalescentProbabilityIndex, double[] result) {
        int errCode = BastaJNIWrapper.INSTANCE.accumulateBastaPartialsGrad(instance,operations, operationCount,
                intervalStarts, intervalCount, intervalLengths, populationSizeIndex, coalescentProbabilityIndex, result);
        if (errCode != 0) {
            throw new BeagleException("accumulateBastaPartialsGrad", errCode);
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

    public int getInstance() {
        return this.instance;
    }

}
