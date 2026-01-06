package dr.evomodel.coalescent.basta;

import dr.evomodel.bigfasttree.BigFastTreeIntervals;
import dr.evomodel.coalescent.SingleTreeGriddedNodesTimeline;
import dr.evomodel.coalescent.basta.ProcessOnCoalescentIntervalDelegate.BranchIntervalOperation;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import java.util.List;

public class PiecewiseConstantPopulationSizeModel extends IntervalSpecificPopulationSizeModel {
    
    private Parameter gridPoints;
    private BigFastTreeIntervals treeIntervals;

    private SingleTreeGriddedNodesTimeline cachedTimeline;
    
    public PiecewiseConstantPopulationSizeModel(String name, 
                                               Parameter populationSizeParameter,
                                               int stateCount, 
                                               int numIntervals) {
        this(name, populationSizeParameter, null, stateCount, numIntervals);
    }
    
    public PiecewiseConstantPopulationSizeModel(String name, 
                                               Parameter populationSizeParameter,
                                               Parameter gridPoints,
                                               int stateCount, 
                                               int numIntervals) {
        super(name, populationSizeParameter, stateCount, numIntervals);
        
        this.gridPoints = gridPoints;
        
        if (gridPoints != null) {
            int numGridSegments = gridPoints.getDimension() + 1;
            int expectedDimWithGrid = stateCount * numGridSegments;
            int paramDim = populationSizeParameter.getDimension();
            
            if (paramDim != expectedDimWithGrid) {
                throw new IllegalArgumentException(
                    "Population size parameter dimension with grid points must be " + 
                    expectedDimWithGrid + " (stateCount=" + stateCount + " × (numGridPoints+1)=" + 
                    numGridSegments + "), but got " + paramDim);
            }
        } else {
            int expectedDim = stateCount * Math.max(1, numIntervals);
            int paramDim = populationSizeParameter.getDimension();
            
            if (paramDim != expectedDim) {
                throw new IllegalArgumentException(
                    "Population size parameter dimension for piecewise constant model must be " + 
                    expectedDim + " (stateCount=" + stateCount + " × numIntervals=" + 
                    Math.max(1, numIntervals) + "), but got " + paramDim);
            }
        }
    }
    
    public void setGridPoints(Parameter gridPoints) {
        this.gridPoints = gridPoints;
        if (gridPoints != null) {
            addVariable(gridPoints);
        }
    }
    
    public void setTreeIntervals(BigFastTreeIntervals treeIntervals) {
        this.treeIntervals = treeIntervals;
    }
    
    public Parameter getGridPoints() {
        return gridPoints;
    }
    
    @Override
    public int getRequiredPopulationSizeStorageSize() {
        if (gridPoints != null) {
            int numGridSegments = gridPoints.getDimension() + 1;
            int gridBasedSize = stateCount + numGridSegments * stateCount;
            int intervalBasedSize = stateCount + numIntervals * stateCount;
            return Math.max(gridBasedSize, intervalBasedSize);
        } else {
            return super.getRequiredPopulationSizeStorageSize();
        }
    }
    
    @Override
    public void extractCombinedSizesAndIntegrals(BastaInternalStorage storage, double[] combinedArray) {
        if (gridPoints != null) {
            int totalSize = getRequiredPopulationSizeStorageSize();

            for (int i = 0; i < totalSize; i++) {
                combinedArray[i] = getCurrentSize(storage, i);
            }

            for (int i = 0; i < totalSize; i++) {
                combinedArray[totalSize + i] = getCurrentIntegral(storage, i);
            }
        } else {
            super.extractCombinedSizesAndIntegrals(storage, combinedArray);
        }
    }
    
    @Override
    protected double calculatePopulationSizeAtTime(int state, int interval,
                                                  double time,
                                                  double intervalStartTime,
                                                  double intervalEndTime) {
        if (gridPoints == null) {
            int paramIndex = interval * stateCount + state;
            return populationSizeParameter.getParameterValue(paramIndex);
        } else {
            int gridIdx = findGridSegmentIndex(gridPoints, time);
            int paramIndex = gridIdx * stateCount + state;
            return populationSizeParameter.getParameterValue(paramIndex);
        }
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
    public PopulationStatistics calculatePopulationStatistics(
            List<Integer> intervalStarts,
            List<BranchIntervalOperation> branchIntervalOperations,
            int stateCount) {
        
        if (gridPoints == null || treeIntervals == null) {
            return super.calculatePopulationStatistics(intervalStarts, branchIntervalOperations, stateCount);
        }
        
        cachedTimeline = new SingleTreeGriddedNodesTimeline(treeIntervals, gridPoints);
        double[] mergedTimeLine = cachedTimeline.getMergedTimeLine();
        int[] gridIndices = cachedTimeline.getGridIndices();
        double[] treeEventTimes = cachedTimeline.getTreeEventTimes();
        
        int numGridSegments = gridPoints.getDimension() + 1;
        int numIntervals = intervalStarts.size() - 1;

        int gridBasedSize = stateCount + numGridSegments * stateCount;
        int integralBasedSize = stateCount + numIntervals * stateCount;
        int requiredStorageSize = Math.max(gridBasedSize, integralBasedSize);
        
        double[] sizes = new double[requiredStorageSize];
        double[] integrals = new double[requiredStorageSize];
        int[] populationSizeIndices = new int[numIntervals];
        int[] integralIndices = new int[numIntervals];

        for (int k = 0; k < stateCount; ++k) {
            sizes[k] = populationSizeParameter.getParameterValue(k);
        }
        
        int currentGridSegment = 0;
        int mergedIdx = 1;
        
        int firstNonZeroIdx = 0;
        while (firstNonZeroIdx < treeEventTimes.length && treeEventTimes[firstNonZeroIdx] == 0.0) {
            firstNonZeroIdx++;
        }

        for (int interval = 0; interval < numIntervals; interval++) {
            double intervalEndTime = treeEventTimes[firstNonZeroIdx + interval];
            double[] integralAccumulator = new double[stateCount];

            while (mergedIdx < mergedTimeLine.length && mergedTimeLine[mergedIdx] <= intervalEndTime) {
                double segmentStart = mergedTimeLine[mergedIdx - 1];
                double segmentEnd = mergedTimeLine[mergedIdx];
                double segmentLength = segmentEnd - segmentStart;

                while (currentGridSegment < gridIndices.length && mergedIdx > gridIndices[currentGridSegment]) {
                    currentGridSegment++;
                }
                
                if (segmentLength > 0) {
                    for (int k = 0; k < stateCount; ++k) {
                        int paramIndex = currentGridSegment * stateCount + k;
                        double popSize = populationSizeParameter.getParameterValue(paramIndex);
                        integralAccumulator[k] += segmentLength / popSize;
                    }
                }
                mergedIdx++;
            }
            
            int intervalOffset = stateCount + interval * stateCount;
            for (int k = 0; k < stateCount; ++k) {
                integrals[intervalOffset + k] = integralAccumulator[k];
                int paramIndex = currentGridSegment * stateCount + k;
                sizes[intervalOffset + k] = populationSizeParameter.getParameterValue(paramIndex);
            }
            
            populationSizeIndices[interval] = intervalOffset;
            integralIndices[interval] = intervalOffset;
        }
        
        return new PopulationStatistics(sizes, integrals, populationSizeIndices, integralIndices, requiredStorageSize);
    }
    
    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        super.handleVariableChangedEvent(variable, index, type);
        cachedTimeline = null;
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
