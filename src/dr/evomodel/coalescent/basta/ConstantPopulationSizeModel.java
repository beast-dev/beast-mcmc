package dr.evomodel.coalescent.basta;

import dr.inference.model.Parameter;
import dr.evomodel.coalescent.basta.ProcessOnCoalescentIntervalDelegate.BranchIntervalOperation;
import java.util.List;


public class ConstantPopulationSizeModel extends AbstractPopulationSizeModel {
    
    private final Parameter populationSizeParameter;
    
    public ConstantPopulationSizeModel(String name, Parameter populationSizeParameter, 
                                     int stateCount, int numIntervals) {
        super(name, stateCount, numIntervals);
        
        this.populationSizeParameter = populationSizeParameter;
        addVariable(populationSizeParameter);
        
        if (populationSizeParameter.getDimension() != stateCount) {
            throw new IllegalArgumentException("Population size parameter dimension must equal state count");
        }
    }
    
    @Override
    protected int calculateArrayLength() {
        return stateCount;
    }
    
    @Override
    public int getPopulationSizeIndex(int interval, int state) {
        return state;
    }
    
    @Override
    public void updatePopulationSizes() {
        for (int i = 0; i < stateCount; i++) {
            populationSizes[i] = populationSizeParameter.getParameterValue(i);
        }
    }
    
    @Override
    public double getIntegralMultiplier(double intervalLength) {
        return intervalLength;
    }
    
    @Override
    public double calculateIntegral(int state, int interval, double intervalStartTime, double intervalLength) {
        if (!populationSizesKnown) {
            updatePopulationSizes();
            populationSizesKnown = true;
        }
        double popSize = populationSizes[state];
        return intervalLength / popSize;
    }
    
    @Override
    public double getPopulationSizeAtTime(int state, int interval, double time) {
        if (!populationSizesKnown) {
            updatePopulationSizes();
            populationSizesKnown = true;
        }
        return populationSizes[state];
    }
    
    @Override
    public PopulationSizeModelType getModelType() {
        return PopulationSizeModelType.CONSTANT;
    }
    
    @Override
    public boolean requiresIntervalSpecificStorage() {
        return false;
    }
    
    public Parameter getPopulationSizeParameter() {
        return populationSizeParameter;
    }
    
    @Override
    public void precalculatePopulationSizesAndIntegrals(
            List<Integer> intervalStarts,
            List<BranchIntervalOperation> branchIntervalOperations,
            BastaInternalStorage storage,
            int stateCount) {

        //storage.flip();
        for (int k = 0; k < stateCount; ++k) {
            double popSize = populationSizeParameter.getParameterValue(k);
            setCurrentSize(storage, k, popSize);
            setCurrentIntegral(storage, k, 1.0 / popSize);
        }

        for (BranchIntervalOperation operation : branchIntervalOperations) {
            operation.populationSizeIndex = 0;
            operation.integralIndex = 0;
        }
    }
}
