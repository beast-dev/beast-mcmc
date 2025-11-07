package dr.evomodel.coalescent.basta;

import dr.inference.model.*;
import dr.evomodel.coalescent.basta.ProcessOnCoalescentIntervalDelegate.BranchIntervalOperation;
import java.util.List;


public abstract class AbstractPopulationSizeModel extends AbstractModel {
    
    protected final int stateCount;
    protected int numIntervals;
    

    protected double[] storedPopulationSizes;
    protected double[] populationSizes;
    protected boolean populationSizesKnown = false;
    
    public AbstractPopulationSizeModel(String name, int stateCount, int numIntervals) {
        super(name);
        this.stateCount = stateCount;
        this.numIntervals = numIntervals;

        int arrayLength = calculateArrayLength();
        this.populationSizes = new double[arrayLength];
        this.storedPopulationSizes = new double[arrayLength];
    }
    

    protected abstract int calculateArrayLength();
    

    public abstract int getPopulationSizeIndex(int interval, int state);
    

    public abstract void updatePopulationSizes();
    

    public abstract double getIntegralMultiplier(double intervalLength);
    

    public abstract double calculateIntegral(int state, int interval, double intervalStartTime, double intervalLength);


    public abstract double getPopulationSizeAtTime(int state, int interval, double time);
    

    public abstract PopulationSizeModelType getModelType();
    

    public abstract boolean requiresIntervalSpecificStorage();

    /**
     * Returns the required storage size for population sizes.
     * For constant models, this is just stateCount.
     * For interval-specific models, this is stateCount + numIntervals * stateCount.
     */
    public int getRequiredPopulationSizeStorageSize() {
        if (requiresIntervalSpecificStorage()) {
            return stateCount + numIntervals * stateCount;
        } else {
            return stateCount;
        }
    }

    public abstract void precalculatePopulationSizesAndIntegrals(
            List<Integer> intervalStarts,
            List<BranchIntervalOperation> branchIntervalOperations,
            BastaInternalStorage storage,
            int stateCount);
    

    public double[] getPopulationSizes() {
        if (!populationSizesKnown) {
            updatePopulationSizes();
            populationSizesKnown = true;
        }
        return populationSizes;
    }
    

    protected void markDirty() {
        populationSizesKnown = false;
        fireModelChanged();
    }
    
    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        markDirty();
    }
    
    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        markDirty();
    }
    
    @Override
    protected void storeState() {

    }
    
    @Override
    protected void restoreState() {

    }
    
    @Override
    protected void acceptState() {

    }


    public enum PopulationSizeModelType {
        CONSTANT(0),
        EXPONENTIAL_GROWTH(1),
        PIECEWISE_CONSTANT(2);

        private final int typeIndex;

        PopulationSizeModelType(int typeIndex) {
            this.typeIndex = typeIndex;
        }

        public int getTypeIndex() {
            return typeIndex;
        }
    }
    

    public void setIntervalCount(int intervalCount) {
        this.numIntervals = intervalCount;
        int arrayLength = calculateArrayLength();
        this.populationSizes = new double[arrayLength];
        this.storedPopulationSizes = new double[arrayLength];
    }
    
    public int getNumIntervals() {
        return numIntervals;
    }

    
    protected double getCurrentSize(BastaInternalStorage storage, int index) {
        int bufferSize = storage.sizes.length / 2;
        int offset = storage.sizesBufferHelper.getOffsetIndex(0) * bufferSize;
        return storage.sizes[offset + index];
    }
    
    protected void setCurrentSize(BastaInternalStorage storage, int index, double value) {
        int bufferSize = storage.sizes.length / 2;
        int offset = storage.sizesBufferHelper.getOffsetIndex(0) * bufferSize;
        storage.sizes[offset + index] = value;
    }
    
    protected double getCurrentIntegral(BastaInternalStorage storage, int index) {
        int bufferSize = storage.integrals.length / 2;
        int offset = storage.integralsBufferHelper.getOffsetIndex(0) * bufferSize;
        return storage.integrals[offset + index];
    }
    
    protected void setCurrentIntegral(BastaInternalStorage storage, int index, double value) {
        int bufferSize = storage.integrals.length / 2;
        int offset = storage.integralsBufferHelper.getOffsetIndex(0) * bufferSize;
        storage.integrals[offset + index] = value;
    }
    
    protected int getCurrentSizesOffset(BastaInternalStorage storage) {
        int bufferSize = storage.sizes.length / 2;
        return storage.sizesBufferHelper.getOffsetIndex(0) * bufferSize;
    }
    
    protected int getCurrentIntegralsOffset(BastaInternalStorage storage) {
        int bufferSize = storage.integrals.length / 2;
        return storage.integralsBufferHelper.getOffsetIndex(0) * bufferSize;
    }
    
    public void extractCombinedSizesAndIntegrals(BastaInternalStorage storage, double[] combinedArray) {
        int requiredStorageSize = getRequiredPopulationSizeStorageSize();

        for (int i = 0; i < requiredStorageSize; i++) {
            combinedArray[i] = getCurrentSize(storage, i);
        }

        for (int i = 0; i < requiredStorageSize; i++) {
            combinedArray[requiredStorageSize + i] = getCurrentIntegral(storage, i);
        }
    }
}
