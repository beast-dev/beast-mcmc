package dr.evomodel.bigfasttree;

import dr.evolution.coalescent.IntervalList;
import dr.evolution.coalescent.IntervalType;
import dr.evolution.util.Units;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

import java.util.Arrays;
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


    public int getSampleEvents(int i) {
        if (i < 0) {
            return sampleCountAt0;
        } else {
            return sampleCounts[i];
        }
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
        double[] tmpWaitTimes = new double[treeIntervals.getIntervalCount()];
        double[] tmpStartTimes = new double[treeIntervals.getIntervalCount()];
        int[] tmpLineageCounts = new int[treeIntervals.getIntervalCount()];
        int[] tmpCoalescentCounts = new int[treeIntervals.getIntervalCount()];
        int[] tmpSampleCounts = new int[treeIntervals.getIntervalCount()];

        // TODO probably off by 1 in modelIntervals


        // TODO what Eldritch horror do we risk summoning if we use this approach?
        int nLineages = 0;
        int idx = 0;
        double treeIntervalStartTime = treeIntervals.getIntervalTime(idx);
        double waitTime = treeIntervals.getInterval(idx);
        double eventTime = treeIntervalStartTime + waitTime;
        // > 1 sample taken at time 0
        if ( eventTime < epsilon ) {
            while (eventTime < epsilon) {
                nLineages++;
                treeIntervalStartTime = treeIntervals.getIntervalTime(idx);
                waitTime = treeIntervals.getInterval(idx);
                eventTime = treeIntervalStartTime + waitTime;
                idx++;
            }
            sampleCountAt0 = nLineages;
            idx--;
        } else {
            nLineages = 1;
            sampleCountAt0 = 1;
            idx = 0;
        }

//        System.err.println("Starting with " + sampleCountAt0 + " lineages");

        double modelIntervalTime = modelIntervals[0];
        tmpStartTimes[0] = treeIntervals.getStartTime();
        int modelIndex = 0;
        int compressedIndex = -1;
        for (int treeIndex = idx; treeIndex < treeIntervals.getIntervalCount(); treeIndex++) {
            treeIntervalStartTime = treeIntervals.getIntervalTime(treeIndex);
            waitTime = treeIntervals.getInterval(treeIndex);
            eventTime = treeIntervalStartTime + waitTime;
//            System.err.println("interval starting at " + treeIntervalStartTime + " with duration " + waitTime + " starting with " + nLineages + " lineages and ending with a " + treeIntervals.getIntervalType(treeIndex));

            if (Math.abs(eventTime - modelIntervalTime) < epsilon) {
                eventTime = modelIntervalTime;
                if (Math.abs(treeIntervalStartTime - modelIntervalTime) < epsilon) {
                    waitTime = 0.0;
                }
            }

            if (eventTime > modelIntervalTime) {
                modelIndex++;
                modelIntervalTime = modelIntervals[modelIndex];
            }

            if (waitTime > 0.0) {
                compressedIndex++;
//                System.err.println("  waitTime > 0.0");
                tmpStartTimes[compressedIndex] = treeIntervalStartTime;
                tmpWaitTimes[compressedIndex] = waitTime;
            }
            
            if (treeIntervals.getIntervalType(treeIndex) == IntervalType.SAMPLE) {
                nLineages++;
                tmpSampleCounts[compressedIndex]++;
//                System.err.println("  tmpSampleCounts[" + compressedIndex + "] = " + tmpSampleCounts[compressedIndex]);
            } else if (treeIntervals.getIntervalType(treeIndex) == IntervalType.COALESCENT) {
                tmpCoalescentCounts[compressedIndex]++;
                nLineages--;
//                System.err.println("  tmpCoalescentCounts[" + compressedIndex + "] = " + tmpCoalescentCounts[compressedIndex]);
            } else {
//                throw new RuntimeException("Tree includes unexpected event type.");
            }
            
            tmpLineageCounts[compressedIndex] = nLineages + tmpCoalescentCounts[compressedIndex] - tmpSampleCounts[compressedIndex];

        }

        // clean up
        intervalCount = compressedIndex + 1;
        intervalsKnown = true;


        startTimes = Arrays.copyOf(tmpStartTimes,intervalCount);
        waitTimes = Arrays.copyOf(tmpWaitTimes,intervalCount);
        lineageCounts = Arrays.copyOf(tmpLineageCounts,intervalCount);
        coalescentCounts = Arrays.copyOf(tmpCoalescentCounts,intervalCount);
        sampleCounts = Arrays.copyOf(tmpSampleCounts,intervalCount);

//        System.err.println("Reporting on times:");
//        System.err.println("start " + new dr.math.matrixAlgebra.Vector(startTimes));
//        System.err.println("waits " + new dr.math.matrixAlgebra.Vector(waitTimes));
//        System.err.println("samps " + new dr.math.matrixAlgebra.Vector(sampleCounts));
//        System.err.println("coals " + new dr.math.matrixAlgebra.Vector(coalescentCounts));
//        System.err.println("linea " + new dr.math.matrixAlgebra.Vector(lineageCounts));
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
    private int sampleCountAt0;
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
