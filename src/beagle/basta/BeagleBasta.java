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

    void getBastaBuffer(int index, double[] buffer);
}
