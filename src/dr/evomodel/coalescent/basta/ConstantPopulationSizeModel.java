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
        return false;
    }
    
    public Parameter getPopulationSizeParameter() {
        return populationSizeParameter;
    }
    
    @Override
    public PopulationStatistics calculatePopulationStatistics(
            List<Integer> intervalStarts,
            List<BranchIntervalOperation> branchIntervalOperations,
            int stateCount) {
        
        int numIntervals = intervalStarts.size() - 1;

        double[] sizes = new double[stateCount];
        double[] integrals = new double[stateCount];
        int[] populationSizeIndices = new int[numIntervals];

        for (int k = 0; k < stateCount; ++k) {
            double popSize = populationSizeParameter.getParameterValue(k);
            sizes[k] = popSize;
            integrals[k] = 1.0 / popSize;
        }

        for (int interval = 0; interval < numIntervals; interval++) {
            populationSizeIndices[interval] = 0;  // All intervals use the same constant population size at index 0
        }
        
        return new PopulationStatistics(sizes, integrals, populationSizeIndices, stateCount);
    }
}
