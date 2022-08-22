package dr.evomodel.bigfasttree;

import dr.evolution.coalescent.IntervalList;
import dr.evolution.coalescent.IntervalType;
import dr.evolution.util.Units;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

import java.util.List;

public class ModelCompressedBigFastTreeIntervals extends AbstractModel implements Units, IntervalList {

    public ModelCompressedBigFastTreeIntervals(BigFastTreeIntervals treeIntervals, double[] modelIntervals, double tolerance) {
        super("modelCompressedBigFastIntervals");
        this.treeIntervals = treeIntervals;
        this.modelIntervals = modelIntervals;
        this.epsilon = tolerance;
        int maxEventCount = treeIntervals.getIntervalCount();

        addModel(treeIntervals);

        calculateIntervals();
    }

    @Override
    public int getIntervalCount() {
        if (!intervalsKnown) {
            calculateIntervals();
        }
        return intervalCount;
    }

    @Override
    public int getSampleCount() {
        return treeIntervals.getSampleCount();
    }

    @Override
    public double getStartTime() {
        return treeIntervals.getStartTime();
    }

    @Override
    public double getInterval(int i) {
        if (!intervalsKnown) {
            calculateIntervals();
        }
        return waitTimes[i];
    }

    @Override
    public double getIntervalTime(int i) {
        if (!intervalsKnown) {
            calculateIntervals();
        }
        return startTimes[i];
    }

    @Override
    public int getLineageCount(int i) {
        if (!intervalsKnown) {
            calculateIntervals();
        }
        return lineageCounts[i];
    }

    @Override
    public int getCoalescentEvents(int i) {
        return coalescentCounts[i];
    }

    @Override
    public IntervalType getIntervalType(int i) {
        if (!intervalsKnown) {
            calculateIntervals();
        }
        return types[i];
    }

    @Override
    public double getTotalDuration() {
        return treeIntervals.getTotalDuration();
    }

    @Override
    public boolean isBinaryCoalescent() {
        return false;
    }

    @Override
    public boolean isCoalescentOnly() {
        return false;
    }

    @Override
    public void calculateIntervals() {
        waitTimes = new double[treeIntervals.getIntervalCount()];
        startTimes = new double[treeIntervals.getIntervalCount()];
        lineageCounts = new int[treeIntervals.getIntervalCount()];
        coalescentCounts = new int[treeIntervals.getIntervalCount()];
        sampleCounts = new int[treeIntervals.getIntervalCount()];

        // TODO probably off by 1 in modelIntervals

        int nLineages = 1;
        double modelIntervalTime = modelIntervals[0];
        startTimes[0] = treeIntervals.getStartTime();

        double treeIntervalStartTime;
        double waitTime;
        double eventTime;


        // handle rest
        int modelIndex = 0;
        int compressedIndex = 0;
        for (int treeIndex = 0; treeIndex < treeIntervals.getIntervalCount(); treeIndex++) {
            treeIntervalStartTime = treeIntervals.getIntervalTime(treeIndex);
            waitTime = treeIntervals.getInterval(treeIndex);
            eventTime = treeIntervalStartTime + waitTime;

            if (eventTime > modelIntervalTime) {
                modelIndex++;
                modelIntervalTime = modelIntervals[modelIndex];
            }

            if (Math.abs(eventTime - modelIntervalTime) < epsilon) {
                eventTime = modelIntervalTime;
            }
            
            if (eventTime > startTimes[compressedIndex]) {
                startTimes[compressedIndex] = treeIntervalStartTime;
                if (compressedIndex == 0) {
                    waitTimes[compressedIndex] = eventTime - treeIntervalStartTime;
                } else {
                    waitTimes[compressedIndex] = treeIntervalStartTime - startTimes[compressedIndex - 1];
                }
            }
            
            if (treeIntervals.getIntervalType(treeIndex) == IntervalType.SAMPLE) {
                nLineages++;
                sampleCounts[compressedIndex]++;
            } else if (treeIntervals.getIntervalType(treeIndex) == IntervalType.COALESCENT) {
                coalescentCounts[compressedIndex]++;
                nLineages--;
            } else {
                throw new RuntimeException("Tree includes unexpected event type.");
            }
            
            lineageCounts[compressedIndex] = nLineages;

            if (eventTime > startTimes[compressedIndex]) {
                compressedIndex++;
            }

//            if (eventTime == startTimes[compressedIndex - 1]) { // this event is concurrent with the last one
//                if (treeIntervals.getIntervalType(treeIndex) == IntervalType.SAMPLE) {
//                    sampleCounts[compressedIndex]++;
//                } else if (treeIntervals.getIntervalType(treeIndex) == IntervalType.COALESCENT) {
//                    coalescentCounts[compressedIndex]++;
//                } else {
//                    throw new RuntimeException("Tree includes unexpected event type.");
//                }
//            } else { // this is a new event
//                startTimes[compressedIndex] = treeIntervalStartTime;
//                waitTimes[compressedIndex] = treeIntervalStartTime - startTimes[compressedIndex - 1];
//                
//                compressedIndex++;
//            }
        }

        // clean up
        double[] tmpIntervals = new double[compressedIndex - 1];
        System.arraycopy(startTimes,0,tmpIntervals,0,compressedIndex - 1);
        startTimes = tmpIntervals;

        System.arraycopy(waitTimes,0,tmpIntervals,0,compressedIndex - 1);
        waitTimes = tmpIntervals;

        int[] tmpCounts = new int[compressedIndex - 1];
        System.arraycopy(sampleCounts,0,tmpCounts,0,compressedIndex - 1);
        sampleCounts = tmpCounts;

        System.arraycopy(coalescentCounts,0,tmpCounts,0,compressedIndex - 1);
        coalescentCounts = tmpCounts;

        System.arraycopy(lineageCounts,0,tmpCounts,0,compressedIndex - 1);
        lineageCounts = tmpCounts;

    }

    private Type units = Type.GENERATIONS;

    @Override
    public final Type getUnits() {
        return units;
    }

    @Override
    public final void setUnits(Type units) {
        this.units = units;
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {

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

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {

    }

    private double[] startTimes;
    private double[] waitTimes;
    private int[] coalescentCounts;
    private int[] sampleCounts;
    private int[] lineageCounts;
    private IntervalType[] types;

    private boolean intervalsKnown;
    private boolean storedIntervalsKnown;

    private boolean onlyUpdateTimes;
    private boolean storedOnlyUpdateTimes;

    private boolean dirty;
    private int intervalCount = 0;

    private BigFastTreeIntervals treeIntervals;
    private double[] modelIntervals;
    private double epsilon;
}
