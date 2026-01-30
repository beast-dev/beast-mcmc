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
    }
    
    @Override
    public boolean requiresIntervalSpecificStorage() {
        return true;
    }

    protected abstract double calculatePopulationSizeAtTime(int state, int interval, 
                                                           double time, 
                                                           double intervalStartTime, 
                                                           double intervalEndTime);

    protected abstract double calculateIntervalIntegral(int state, int interval,
                                                       double intervalStartTime,
                                                       double intervalLength);

    protected void getBaseSizes(double[] baseSizes) {
        for (int k = 0; k < stateCount; ++k) {
            baseSizes[k] = populationSizeParameter.getParameterValue(k);
        }
    }

    @Override
    public PopulationStatistics calculatePopulationStatistics(
            List<Integer> intervalStarts,
            List<BranchIntervalOperation> branchIntervalOperations,
            double[] intervalLengths,
            int stateCount) {
        
        int numIntervals = intervalStarts.size() - 1;
        int requiredStorageSize = stateCount + numIntervals * stateCount;
        
        double[] sizes = new double[requiredStorageSize];
        double[] integrals = new double[requiredStorageSize];
        int[] populationSizeIndices = new int[numIntervals];

        getBaseSizes(sizes);
        
        double intervalStartTime = 0.0;
        for (int interval = 0; interval < numIntervals; ++interval) {
            double intervalLength = intervalLengths[interval];
            double intervalEndTime = intervalStartTime + intervalLength;
            
            int populationSizeIndex = stateCount + interval * stateCount;
            
            for (int k = 0; k < stateCount; ++k) {
                double popSizeAtEnd = calculatePopulationSizeAtTime(k, interval, 
                    intervalEndTime, intervalStartTime, intervalEndTime);
                sizes[populationSizeIndex + k] = popSizeAtEnd;
                
                double integral = calculateIntervalIntegral(k, interval, intervalStartTime, intervalLength);
                integrals[populationSizeIndex + k] = integral;
            }
            
            populationSizeIndices[interval] = populationSizeIndex;
            
            intervalStartTime = intervalEndTime;
        }
        
        return new PopulationStatistics(sizes, integrals, populationSizeIndices, requiredStorageSize);
    }

    public Parameter getPopulationSizeParameter() {
        return populationSizeParameter;
    }
}
