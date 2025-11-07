package dr.evomodel.coalescent.basta;

import dr.inference.model.Parameter;
import dr.evomodel.coalescent.basta.ProcessOnCoalescentIntervalDelegate.BranchIntervalOperation;
import java.util.List;

public abstract class IntervalSpecificPopulationSizeModel extends AbstractPopulationSizeModel {
    
    protected final Parameter populationSizeParameter;
    
    public IntervalSpecificPopulationSizeModel(String name, Parameter populationSizeParameter,
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
        return stateCount + numIntervals * stateCount;
    }
    
    @Override
    public boolean requiresIntervalSpecificStorage() {
        return true;
    }
    
    @Override
    public int getPopulationSizeIndex(int interval, int state) {
        return stateCount + interval * stateCount + state;
    }
    

    protected abstract double calculatePopulationSizeAtTime(int state, int interval, 
                                                           double time, 
                                                           double intervalStartTime, 
                                                           double intervalEndTime);
    

    protected abstract double calculateIntervalIntegral(int state, int interval,
                                                       double intervalStartTime,
                                                       double intervalLength);


    protected abstract void storeBaseSizes(BastaInternalStorage storage);
    
    @Override
    public final void precalculatePopulationSizesAndIntegrals(
            List<Integer> intervalStarts,
            List<BranchIntervalOperation> branchIntervalOperations,
            BastaInternalStorage storage,
            int stateCount) {
        
        int numIntervals = intervalStarts.size() - 1;

        storeBaseSizes(storage);

        double intervalStartTime = 0.0;
        for (int interval = 0; interval < numIntervals; ++interval) {
            int start = intervalStarts.get(interval);
            int end = intervalStarts.get(interval + 1);
            
            double intervalLength = branchIntervalOperations.get(start).intervalLength;
            double intervalEndTime = intervalStartTime + intervalLength;
            
            int intervalOffset = interval * stateCount;

            for (int k = 0; k < stateCount; ++k) {
                double popSizeAtEnd = calculatePopulationSizeAtTime(k, interval, 
                    intervalEndTime, intervalStartTime, intervalEndTime);
                setCurrentSize(storage, stateCount + intervalOffset + k, popSizeAtEnd);

                double integral = calculateIntervalIntegral(k, interval, 
                    intervalStartTime, intervalLength);
                int integralOffset = stateCount + intervalOffset;
                setCurrentIntegral(storage, integralOffset + k, integral);
            }

            int populationSizeIndex = stateCount + interval * stateCount;
            for (int i = start; i < end; ++i) {
                branchIntervalOperations.get(i).populationSizeIndex = populationSizeIndex;
            }
            
            intervalStartTime = intervalEndTime;
        }
    }
    
    public Parameter getPopulationSizeParameter() {
        return populationSizeParameter;
    }
}

