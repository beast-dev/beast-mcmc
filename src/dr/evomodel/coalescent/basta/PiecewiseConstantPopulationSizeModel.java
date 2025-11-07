package dr.evomodel.coalescent.basta;

import dr.inference.model.Parameter;

public class PiecewiseConstantPopulationSizeModel extends IntervalSpecificPopulationSizeModel {
    
    public PiecewiseConstantPopulationSizeModel(String name, 
                                               Parameter populationSizeParameter,
                                               int stateCount, 
                                               int numIntervals) {
        super(name, populationSizeParameter, stateCount, numIntervals);

        int expectedDim = stateCount * Math.max(1, numIntervals);
        int paramDim = populationSizeParameter.getDimension();
        
        if (paramDim != expectedDim) {
            throw new IllegalArgumentException(
                "Population size parameter dimension for piecewise constant model must be " + 
                expectedDim + " (stateCount=" + stateCount + " × numIntervals=" + 
                Math.max(1, numIntervals) + "), but got " + paramDim);
        }
    }
    
    @Override
    protected void storeBaseSizes(BastaInternalStorage storage) {
        for (int k = 0; k < stateCount; ++k) {
            setCurrentSize(storage, k, populationSizeParameter.getParameterValue(k));
        }
    }
    
    @Override
    protected double calculatePopulationSizeAtTime(int state, int interval,
                                                  double time,
                                                  double intervalStartTime,
                                                  double intervalEndTime) {
        int paramIndex = interval * stateCount + state;
        return populationSizeParameter.getParameterValue(paramIndex);
    }
    
    @Override
    protected double calculateIntervalIntegral(int state, int interval,
                                              double intervalStartTime,
                                              double intervalLength) {
        int paramIndex = interval * stateCount + state;
        double size = populationSizeParameter.getParameterValue(paramIndex);
        return intervalLength / size;
    }
    
    @Override
    public void updatePopulationSizes() {
        for (int interval = 0; interval < numIntervals; interval++) {
            for (int state = 0; state < stateCount; state++) {
                int paramIndex = interval * stateCount + state;
                int storageIndex = getPopulationSizeIndex(interval, state);
                populationSizes[storageIndex] = populationSizeParameter.getParameterValue(paramIndex);
            }
        }

        for (int state = 0; state < stateCount; state++) {
            populationSizes[state] = populationSizeParameter.getParameterValue(state);
        }
    }
    
    @Override
    public double getIntegralMultiplier(double intervalLength) {
        return 1.0;
    }
    
    @Override
    public double calculateIntegral(int state, int interval, double intervalStartTime, double intervalLength) {
        return calculateIntervalIntegral(state, interval, intervalStartTime, intervalLength);
    }
    
    @Override
    public double getPopulationSizeAtTime(int state, int interval, double time) {
        return calculatePopulationSizeAtTime(state, interval, time, 0, time);
    }
    
    @Override
    public PopulationSizeModelType getModelType() {
        return PopulationSizeModelType.PIECEWISE_CONSTANT;
    }
}
