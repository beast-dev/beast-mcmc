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

    void updateBastaPartialsGrad(int[] operations, int operationCount, int[] intervals, int intervalCount, int populationSizeIndex, int coalescentProbabilityIndex);

    void updateTransitionMatricesGrad(int[] transitionMatrixIndices, double[] branchLengths, int count);

    void accumulateBastaPartialsGrad(int[] operations, int operationCount,  final int[] intervalStarts, int intervalCount, final double[] intervalLengths, int populationSizeIndex, int coalescentProbabilityIndex, double[] result);
}
