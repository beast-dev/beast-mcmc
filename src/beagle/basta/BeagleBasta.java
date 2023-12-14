package beagle.basta;

import beagle.Beagle;

public interface BeagleBasta extends Beagle {

    void updateBastaPartials(final int[] operations,
                             int operationCount,
                             int populationSizeIndex);

    void accumulateBastaPartials(final int[] operations,
                                 int operationCount,
                                 final int[] segments,
                                 int segmentCount);
}
