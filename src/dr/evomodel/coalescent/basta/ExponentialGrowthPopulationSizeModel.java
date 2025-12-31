package dr.evomodel.coalescent.basta;

import dr.inference.model.Parameter;

public class ExponentialGrowthPopulationSizeModel extends IntervalSpecificPopulationSizeModel {
    
    private final Parameter growthRateParameter;
    
    public ExponentialGrowthPopulationSizeModel(String name, 
                                              Parameter populationSizeParameter,
                                              Parameter growthRateParameter,
                                              int stateCount, 
                                              int numIntervals) {
        super(name, populationSizeParameter, stateCount, numIntervals);
        
        this.growthRateParameter = growthRateParameter;
        addVariable(growthRateParameter);

        int popSizesDim = populationSizeParameter.getDimension();
        if (!(popSizesDim == 1 || popSizesDim == stateCount)) {
            throw new IllegalArgumentException("Population size parameter dimension must equal state count for exponential growth model, got " +
                populationSizeParameter.getDimension() + " but expected " + stateCount);
        }

        int growthRateDim = growthRateParameter.getDimension();
        if (!(growthRateDim == 1 || growthRateDim == stateCount)) {
            throw new IllegalArgumentException("Growth rate parameter dimension must equal state count");
        }
    }
    
    @Override
    protected void storeBaseSizes(BastaInternalStorage storage) {
        for (int k = 0; k < stateCount; ++k) {
            int j = populationSizeParameter.getDimension() == 1 ? 0 : k;
            setCurrentSize(storage, k, populationSizeParameter.getParameterValue(j));
            storage.rates[k] = growthRateParameter.getParameterValue(j);
        }
    }
    
    @Override
    protected double calculatePopulationSizeAtTime(int state, int interval,
                                                  double time,
                                                  double intervalStartTime,
                                                  double intervalEndTime) {
        int j = populationSizeParameter.getDimension() == 1 ? 0 : state;
        double N0 = populationSizeParameter.getParameterValue(j);
        double r = growthRateParameter.getParameterValue(j);
        
        if (N0 <= 0.0) {
            throw new RuntimeException("Population size must be positive, got N0=" + N0);
        }
        
        // N(t) = N0 * exp(-r*t)
        return N0 * Math.exp(-r * time);
    }
    
    @Override
    protected double calculateIntervalIntegral(int state, int interval,
                                              double intervalStartTime,
                                              double intervalLength) {
        int j = populationSizeParameter.getDimension() == 1 ? 0 : state;
        double N0 = populationSizeParameter.getParameterValue(j);
        double r = growthRateParameter.getParameterValue(j);
        double intervalEndTime = intervalStartTime + intervalLength;
        
        if (Math.abs(r) > 1e-10) {
            double maxExp = Math.max(r * intervalStartTime, r * intervalEndTime);
            if (maxExp > 700.0) {
                throw new RuntimeException(String.format(
                    "Exponential overflow in integral calculation: r=%g, time=%g, maxExp=%g",
                    r, intervalEndTime, maxExp));
            }
            
            // Integral of 1/N(t) from t1 to t2 = (exp(r*t2) - exp(r*t1)) / (N0 * r)
            double expStartTerm = Math.exp(r * intervalStartTime);
            double expEndTerm = Math.exp(r * intervalEndTime);
            return (expEndTerm - expStartTerm) / (N0 * r);
        } else {
            return intervalLength / N0;
        }
    }
    
    @Override
    public void updatePopulationSizes() {
        for (int i = 0; i < stateCount; i++) {
            int j = populationSizeParameter.getDimension() == 1 ? 0 : i;
            populationSizes[i] = populationSizeParameter.getParameterValue(j);
        }
        

        // HERE

        for (int interval = 0; interval < numIntervals; interval++) {
            for (int state = 0; state < stateCount; state++) {
                int index = getPopulationSizeIndex(interval, state);
                int j = populationSizeParameter.getDimension() == 1 ? 0 : state;
                populationSizes[index] = populationSizeParameter.getParameterValue(j);
            }
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
        return PopulationSizeModelType.EXPONENTIAL_GROWTH;
    }
    
    public Parameter getGrowthRateParameter() {
        return growthRateParameter;
    }
}
