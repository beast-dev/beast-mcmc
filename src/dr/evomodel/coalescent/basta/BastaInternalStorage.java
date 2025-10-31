package dr.evomodel.coalescent.basta;

import dr.evomodel.substmodel.EigenDecomposition;
import dr.evomodel.treedatalikelihood.BufferIndexHelper;

/**
 * @author Yucai Shao
 * @author Marc A. Suchard
 */

public class BastaInternalStorage {

    double[] partials;
    double[] matrices;
    double[] coalescent;

    double[] storedMatrices;

    double[] e;
    double[] f;
    double[] g;
    double[] h;

    double[] sizes;
    double[] integrals;
    final double[] rates;
    
    private final BufferIndexHelper sizesBufferHelper;
    private final BufferIndexHelper integralsBufferHelper;

    final EigenDecomposition[] decompositions; // TODO flatten?

    private int currentNumCoalescentIntervals;
    private int currentNumPartials;

    final private int stateCount;

    public BastaInternalStorage(int maxNumCoalescentIntervals, int treeNodeCount, int stateCount) {

        this.currentNumPartials = 0;
        this.currentNumCoalescentIntervals = 0;

        this.stateCount = stateCount;
        this.sizesBufferHelper = new SimpleBufferIndexHelper();
        this.integralsBufferHelper = new SimpleBufferIndexHelper();

        int baseLength = stateCount;
        this.sizes = new double[baseLength * 2];      // Double for buffer flipping
        this.integrals = new double[baseLength * 2];  // Double for buffer flipping
        this.rates = new double[2 * stateCount];
        this.decompositions = new EigenDecomposition[1];

        resize(3 * treeNodeCount, maxNumCoalescentIntervals, null);
    }

    public void resizeSizesAndIntegrals(int requiredSize) {
        int actualSize = requiredSize * 2;
        if (sizes.length != actualSize) {
            double[] oldSizes = this.sizes;
            double[] oldIntegrals = this.integrals;
            
            this.sizes = new double[actualSize];
            this.integrals = new double[actualSize];

            if (oldSizes != null && oldSizes.length > 0) {
                int copyLength = Math.min(oldSizes.length, actualSize);
                System.arraycopy(oldSizes, 0, this.sizes, 0, copyLength);
                System.arraycopy(oldIntegrals, 0, this.integrals, 0, copyLength);
            }
        }
    }

    static private int getStartingPartialsCount(int maxNumCoalescentIntervals, int treeNodeCount) {
        return maxNumCoalescentIntervals * (treeNodeCount + 1); // TODO much too large
    }

    public void resize(int newNumPartials, int newNumCoalescentIntervals, BastaLikelihood likelihood) {
        resize(newNumPartials, newNumCoalescentIntervals, likelihood, null, 0);
    }
    
    public void resize(int newNumPartials, int newNumCoalescentIntervals, BastaLikelihood likelihood,
                       AbstractPopulationSizeModel populationSizeModel, int numIntervals) {

        if (newNumPartials > currentNumPartials) {
            this.partials = new double[newNumPartials * stateCount];
            this.currentNumPartials =  newNumPartials;
            if (likelihood != null) {
                likelihood.setTipData();
            }
        }

        if (newNumCoalescentIntervals > this.currentNumCoalescentIntervals) {
            this.matrices = new double[newNumCoalescentIntervals * stateCount * stateCount]; // TODO much too small (except for strict-clock)
            this.storedMatrices = new double[newNumCoalescentIntervals * stateCount * stateCount];
            this.coalescent = new double[newNumCoalescentIntervals];

            this.e = new double[newNumCoalescentIntervals * stateCount];
            this.f = new double[newNumCoalescentIntervals * stateCount];
            this.g = new double[newNumCoalescentIntervals * stateCount];
            this.h = new double[newNumCoalescentIntervals * stateCount];

            this.currentNumCoalescentIntervals = newNumCoalescentIntervals;
        }

        if (populationSizeModel != null && numIntervals > 0) {
            if (populationSizeModel.getNumIntervals() != numIntervals) {
                populationSizeModel.setIntervalCount(numIntervals);
            }
            int requiredSize = populationSizeModel.getRequiredPopulationSizeStorageSize();
            resizeSizesAndIntegrals(requiredSize);
        }
    }

    public void storeState() {
        System.arraycopy(matrices, 0, storedMatrices, 0, currentNumCoalescentIntervals * stateCount * stateCount);
        sizesBufferHelper.storeState();
        integralsBufferHelper.storeState();
    }

    public void restoreState() {
        double[] temp = matrices;
        matrices = storedMatrices;
        storedMatrices = temp;
        sizesBufferHelper.restoreState();
        integralsBufferHelper.restoreState();
    }

    public int getCurrentSizesBuffer() {
        int bufferSize = sizes.length / 2;
        int offset = sizesBufferHelper.getOffsetIndex(0) * bufferSize;
        return offset;
    }

    public double getCurrentSize(int index) {
        int bufferSize = sizes.length / 2;
        int offset = sizesBufferHelper.getOffsetIndex(0) * bufferSize;
        return sizes[offset + index];
    }
    
    public void setCurrentSize(int index, double value) {
        int bufferSize = sizes.length / 2;
        int offset = sizesBufferHelper.getOffsetIndex(0) * bufferSize;
        sizes[offset + index] = value;
    }
    
    public double getCurrentIntegral(int index) {
        int bufferSize = integrals.length / 2;
        int offset = integralsBufferHelper.getOffsetIndex(0) * bufferSize;
        return integrals[offset + index];
    }
    
    public void setCurrentIntegral(int index, double value) {
        int bufferSize = integrals.length / 2;
        int offset = integralsBufferHelper.getOffsetIndex(0) * bufferSize;
        if (offset + index >= integrals.length) {
            throw new ArrayIndexOutOfBoundsException("Index " + (offset + index) + 
                " out of bounds for length " + integrals.length + 
                " (buffer offset=" + offset + ", index=" + index + ", bufferSize=" + bufferSize + ")");
        }
        integrals[offset + index] = value;
    }
    
    public int getActualStorageSize() {
        return sizes.length / 2;
    }

    public void flip() {
        sizesBufferHelper.flipOffset(0);
        integralsBufferHelper.flipOffset(0);
    }
    

    public int getCurrentSizesOffset() {
        int bufferSize = sizes.length / 2;
        return sizesBufferHelper.getOffsetIndex(0) * bufferSize;
    }
    

    public double[] getCurrentIntegralsArray() {
        return integrals;
    }

    public int getCurrentIntegralsOffset() {
        int bufferSize = integrals.length / 2;
        return integralsBufferHelper.getOffsetIndex(0) * bufferSize;
    }

    private static class SimpleBufferIndexHelper extends BufferIndexHelper {
        public SimpleBufferIndexHelper() {
            super(1, 0);
        }

        @Override
        protected int computeOffset(int offset) { 
            return offset;
        }
    }
}
