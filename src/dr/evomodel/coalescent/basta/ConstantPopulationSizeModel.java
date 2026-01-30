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
    public PopulationSizeModelType getModelType() {
        return PopulationSizeModelType.CONSTANT;
    }
    
    @Override
    public boolean requiresIntervalSpecificStorage() {
        return true;
    }
    
    public Parameter getPopulationSizeParameter() {
        return populationSizeParameter;
    }
    
    @Override
    public PopulationStatistics calculatePopulationStatistics(
            List<Integer> intervalStarts,
            List<BranchIntervalOperation> branchIntervalOperations,
            double[] intervalLengths,
            int stateCount) {
        
        int numIntervals = intervalStarts.size() - 1;
        int storageSize = stateCount + numIntervals * stateCount;

        double[] sizes = new double[storageSize];
        double[] integrals = new double[storageSize];
        int[] populationSizeIndices = new int[numIntervals];

        for (int interval = 0; interval < numIntervals; interval++) {
            int offset = stateCount + interval * stateCount;
            double length = intervalLengths[interval];
            
            for (int k = 0; k < stateCount; ++k) {
                double popSize = populationSizeParameter.getParameterValue(k);
                sizes[offset + k] = popSize;
                integrals[offset + k] = length / popSize;
            }
            populationSizeIndices[interval] = offset;
        }
        
        return new PopulationStatistics(sizes, integrals, populationSizeIndices, storageSize);
    }
}
