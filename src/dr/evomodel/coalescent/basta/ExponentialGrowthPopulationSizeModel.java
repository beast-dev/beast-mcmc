package dr.evomodel.coalescent.basta;

import dr.inference.model.Parameter;
import dr.evomodel.coalescent.basta.ProcessOnCoalescentIntervalDelegate.BranchIntervalOperation;
import java.util.List;


public class ExponentialGrowthPopulationSizeModel extends AbstractPopulationSizeModel {
    
    private final Parameter populationSizeParameter;
    private final Parameter growthRateParameter;
    
    public ExponentialGrowthPopulationSizeModel(String name, 
                                              Parameter populationSizeParameter,
                                              Parameter growthRateParameter,
                                              int stateCount, 
                                              int numIntervals) {
        super(name, stateCount, numIntervals);
        
        this.populationSizeParameter = populationSizeParameter;
        this.growthRateParameter = growthRateParameter;
        
        addVariable(populationSizeParameter);
        addVariable(growthRateParameter);
        
        if (populationSizeParameter.getDimension() != stateCount) {
            throw new IllegalArgumentException("Population size parameter dimension must equal state count");
        }
        
        if (growthRateParameter.getDimension() != stateCount) {
            throw new IllegalArgumentException("Growth rate parameter dimension must equal state count");
        }

    }
    
    @Override
    protected int calculateArrayLength() {

        return (numIntervals > 0) ? (stateCount + numIntervals * stateCount) : stateCount;
    }
    
    @Override
    public int getPopulationSizeIndex(int interval, int state) {
        return stateCount + interval * stateCount + state;
    }
    
    @Override
    public void updatePopulationSizes() {
        for (int i = 0; i < stateCount; i++) {
            populationSizes[i] = populationSizeParameter.getParameterValue(i);
        }
        

        for (int interval = 0; interval < numIntervals; interval++) {
            for (int state = 0; state < stateCount; state++) {
                int index = getPopulationSizeIndex(interval, state);
                populationSizes[index] = populationSizeParameter.getParameterValue(state);
            }
        }
    }
    
    @Override
    public double getIntegralMultiplier(double intervalLength) {
        return 1.0;
    }
    
    @Override
    public double calculateIntegral(int state, int interval, double intervalStartTime, double intervalLength) {
        if (!populationSizesKnown) {
            updatePopulationSizes();
            populationSizesKnown = true;
        }
        
        double N0 = populationSizeParameter.getParameterValue(state);
        double r = growthRateParameter.getParameterValue(state);
        double intervalEndTime = intervalStartTime + intervalLength;
        
        if (Math.abs(r) < 1e-10) {
            return intervalLength / N0;
        } else {
            // For exponential growth: N(t) = N0 * exp(-r*t)
            // Integral of 1/N(t) from t1 to t2 = (exp(r*t2) - exp(r*t1)) / (N0 * r)
            double expStartTerm = Math.exp(r * intervalStartTime);
            double expEndTerm = Math.exp(r * intervalEndTime);
            return (expEndTerm - expStartTerm) / (N0 * r);
        }
    }
    
    @Override
    public double getPopulationSizeAtTime(int state, int interval, double time) {
        if (!populationSizesKnown) {
            updatePopulationSizes();
            populationSizesKnown = true;
        }
        
        double N0 = populationSizeParameter.getParameterValue(state);
        double r = growthRateParameter.getParameterValue(state);

        return N0 * Math.exp(-r * time);
    }
    
    @Override
    public PopulationSizeModelType getModelType() {
        return PopulationSizeModelType.EXPONENTIAL_GROWTH;
    }
    
    @Override
    public boolean requiresIntervalSpecificStorage() {
        return true;
    }
    
    public Parameter getPopulationSizeParameter() {
        return populationSizeParameter;
    }
    
    public Parameter getGrowthRateParameter() {
        return growthRateParameter;
    }

    
    @Override
    public void precalculatePopulationSizesAndIntegrals(
            List<Integer> intervalStarts,
            List<BranchIntervalOperation> branchIntervalOperations,
             BastaInternalStorage storage,
             int stateCount) {
 
         int numIntervals = intervalStarts.size() - 1;

        // Store base population sizes and rates using helper methods
        for (int k = 0; k < stateCount; ++k) {
            storage.setCurrentSize(k, populationSizeParameter.getParameterValue(k));
            storage.rates[k] = growthRateParameter.getParameterValue(k);
        }
        

        double intervalStartTime = 0.0;
        // Only store interval-specific values for intervals 0 to numIntervals-1
        for (int interval = 0; interval < numIntervals; ++interval) {
            int start = intervalStarts.get(interval);
            int end = intervalStarts.get(interval + 1);

            double intervalLength = branchIntervalOperations.get(start).intervalLength;
            double intervalEndTime = intervalStartTime + intervalLength;
            
            int intervalOffset = interval * stateCount;

            for (int k = 0; k < stateCount; ++k) {
                double N0 = storage.getCurrentSize(k);
                double r = storage.rates[k];

                if (N0 <= 0.0) {
                    throw new RuntimeException("Population size must be positive, got N0=" + N0);
                }

                double popSizeAtEnd = N0 * Math.exp(-r * intervalEndTime);
                storage.setCurrentSize(stateCount + intervalOffset + k, popSizeAtEnd);

                int integralOffset = stateCount + intervalOffset;
                if (Math.abs(r) > 1e-10) {
                    double maxExp = Math.max(r * intervalStartTime, r * intervalEndTime);
                    if (maxExp > 700.0) {
                        throw new RuntimeException(String.format(
                            "Exponential overflow in integral calculation: r=%g, time=%g, maxExp=%g. ",
                            r, intervalEndTime, maxExp));
                    }
                    double expStartTerm = Math.exp(r * intervalStartTime);
                    double expEndTerm = Math.exp(r * intervalEndTime);
                    storage.setCurrentIntegral(integralOffset + k, (expEndTerm - expStartTerm) / (N0 * r));

                } else {
                    storage.setCurrentIntegral(integralOffset + k, intervalLength / N0);
                }
            }

            int populationSizeIndex = stateCount + interval * stateCount;
            for (int i = start; i < end; ++i) {
                branchIntervalOperations.get(i).populationSizeIndex = populationSizeIndex;
            }

            intervalStartTime = intervalEndTime;
        }
    }

}
