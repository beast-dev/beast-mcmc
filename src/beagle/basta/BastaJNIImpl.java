package beagle.basta;

import beagle.BeagleException;
import beagle.BeagleJNIImpl;

public class BastaJNIImpl extends BeagleJNIImpl implements BeagleBasta {

    public BastaJNIImpl(int tipCount,
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

    }

    @Override
    public void updateBastaPartials(int[] operations, int operationCount, int populationSizeIndex) {
        int errCode = BastaJNIWrapper.INSTANCE.updateBastaPartials(instance, operations, operationCount,
                populationSizeIndex);
        if (errCode != 0) {
            throw new BeagleException("updateBastaPartials", errCode);
        }
    }

    @Override
    public void accumulateBastaPartials(int[] operations, int operationCount, int[] segments, int segmentCount) {
        int errCode = BastaJNIWrapper.INSTANCE.accumulateBastaPartials(instance,operations, operationCount,
                segments, segmentCount);
        if (errCode != 0) {
            throw new BeagleException("accumulateBastaPartials", errCode);
        }
    }
}
