package dr.evomodel.coalescent.basta;

import dr.inference.model.*;
import dr.evomodel.coalescent.basta.ProcessOnCoalescentIntervalDelegate.BranchIntervalOperation;
import java.util.List;


public abstract class AbstractPopulationSizeModel extends AbstractModel {
    
    protected final int stateCount;
    protected int numIntervals;
    
    public AbstractPopulationSizeModel(String name, int stateCount, int numIntervals) {
        super(name);
        this.stateCount = stateCount;
        this.numIntervals = numIntervals;
    }

    public static class PopulationStatistics {
        public final double[] sizes;
        public final double[] integrals;
        public final int[] populationSizeIndices;  // Used for indexing both sizes and integrals
        public final int requiredStorageSize;
        
        public PopulationStatistics(double[] sizes, double[] integrals,
                                   int[] populationSizeIndices,
                                   int requiredStorageSize) {
            this.sizes = sizes;
            this.integrals = integrals;
            this.populationSizeIndices = populationSizeIndices;
            this.requiredStorageSize = requiredStorageSize;
        }
    }

    public abstract PopulationStatistics calculatePopulationStatistics(
            List<Integer> intervalStarts,
            List<BranchIntervalOperation> branchIntervalOperations,
            int stateCount);

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
    
    public void setIntervalCount(int intervalCount) {
        this.numIntervals = intervalCount;
    }
    
    public int getNumIntervals() {
        return numIntervals;
    }

    protected void markDirty() {
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
}
