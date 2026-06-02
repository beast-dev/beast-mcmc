/*
 * Intervals.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 *
 */

package dr.evolution.coalescent;

import dr.math.MathUtils;

import java.util.Arrays;

/**
 * A concrete class for a set of coalescent intervals.
 *
 * @author Andrew Rambaut
 */
public class FastIntervals implements MutableIntervalList {


    public FastIntervals(int sampleCount, int coalescentCount) {
        startTime = Double.POSITIVE_INFINITY;
        finishTime = Double.NEGATIVE_INFINITY;

        this.sampleCount = sampleCount;
        this.coalescentCount = coalescentCount;
        eventCount = sampleCount + coalescentCount;
        intervalCount = eventCount - 1;

        sampleTimes = new double[sampleCount];
        coalescentTimes = new double[coalescentCount];
        eventTimes = new double[eventCount];

        intervals = new double[intervalCount];
        intervalTypes = new IntervalType[intervalCount];
        lineageCounts = new int[intervalCount];

        cumulativeSampleCount = 0;
        cumulativeCoalescentCount = 0;

        intervalsKnown = false;
    }
    @Override
    public void copyIntervals(MutableIntervalList intervalList) {
        FastIntervals source = (FastIntervals)intervalList;
        intervalsKnown = source.intervalsKnown;

        assert eventCount == source.eventCount;
        assert sampleCount == source.sampleCount;
        assert coalescentCount == source.coalescentCount;
        assert intervalCount == source.intervalCount;

        startTime = source.startTime;

        if (intervalsKnown) {
            System.arraycopy(source.intervals, 0, intervals, 0, intervals.length);
            System.arraycopy(source.intervalTypes, 0, intervalTypes, 0, intervals.length);
            System.arraycopy(source.lineageCounts, 0, lineageCounts, 0, intervals.length);
        }
    }

    @Override
    public void resetEvents() {
        startTime = Double.POSITIVE_INFINITY;
        finishTime = Double.NEGATIVE_INFINITY;

        intervalsKnown = false;
        cumulativeSampleCount = 0;
        cumulativeCoalescentCount = 0;
    }

    @Override
    public void addSampleEvent(double time) {
        if (time < startTime) {
            startTime = time;
        }

        sampleTimes[cumulativeSampleCount] = time;
        cumulativeSampleCount++;
        intervalsKnown = false;
    }

    @Override
    public void addSampleEvent(double time, int nodeNumber) {
        // ignores the node number
        addSampleEvent(time);
    }

    @Override
    public void addCoalescentEvent(double time) {
        if (time > finishTime) {
            finishTime = time;
        }

        coalescentTimes[cumulativeCoalescentCount] = time;
        cumulativeCoalescentCount++;
        intervalsKnown = false;
    }

    @Override
    public void addCoalescentEvent(double time, int nodeNumber) {
        // ignores the node number
        addCoalescentEvent(time);
    }

    @Override
    public void addMigrationEvent(double time, int destination) {
        throw new UnsupportedOperationException("not supported in FastIntervals");
    }

    @Override
    public void addNothingEvent(double time) {
        throw new UnsupportedOperationException("not supported in FastIntervals");
    }

    @Override
    public int getNodeForEvent(int i) {
        throw new UnsupportedOperationException("not supported in FastIntervals");
    }

    public int getSampleCount() {
        return sampleCount;
    }

    public int getIntervalCount() {
        if (!intervalsKnown) {
            calculateIntervals();
        }
        return intervalCount;
    }

    public double getInterval(int i) {
        if (!intervalsKnown) {
            calculateIntervals();
        }
        return intervals[i];
    }

    public double getIntervalTime(int i){
        if (!intervalsKnown){
            calculateIntervals();
        }
        return eventTimes[i];
    }

    public int getLineageCount(int i) {
        if (!intervalsKnown) {
            calculateIntervals();
        }
        return lineageCounts[i];
    }

    public int getCoalescentEvents(int i) {
        if (!intervalsKnown) {
            calculateIntervals();
        }
        if (i < intervalCount - 1) {
            return lineageCounts[i] - lineageCounts[i + 1];
        } else {
            return lineageCounts[i] - 1;
        }
    }

    public double getStartTime() {
        if (!intervalsKnown) {
            calculateIntervals();
        }
        return startTime;
    }

    public IntervalType getIntervalType(int i) {
        if (!intervalsKnown) {
            calculateIntervals();
        }
        return intervalTypes[i];
    }
    public double getTotalDuration() {
        if (!intervalsKnown) {
            calculateIntervals();
        }
        return finishTime - startTime;
    }

    public boolean isBinaryCoalescent() {
        return true;
    }

    public boolean isCoalescentOnly() {
        return true;
    }

    public void calculateIntervals() {

        if (eventCount < 2) {
            throw new IllegalArgumentException("Too few events to construct intervals");
        }

        // sampleTimes only changes if there is tipdate sampling but it would be complicated
        // to check this. The sort function on an already sorted array is very fast.
        Arrays.sort(sampleTimes);
        Arrays.sort(coalescentTimes);

        // sample should come first
        double lastTime = sampleTimes[0];

        int s = 1;
        int c = 0;
        int i = 0;
        int lineages = 1;

        while (s < sampleCount || c < coalescentCount) {
            if (s < sampleCount && sampleTimes[s] <= coalescentTimes[c]) {
                intervals[i] = sampleTimes[s] - lastTime;
                eventTimes[i] = sampleTimes[s];
                intervalTypes[i] = IntervalType.SAMPLE;
                lineageCounts[i] = lineages;
                lastTime = sampleTimes[s];
                lineages++;
                s++;
            } else {
                intervals[i] = coalescentTimes[c] - lastTime;
                eventTimes[i] = coalescentTimes[c];
                intervalTypes[i] = IntervalType.COALESCENT;
                lineageCounts[i] = lineages;
                lastTime = coalescentTimes[c];
                lineages--;
                c++;
            }
            i++;
        }

        intervalsKnown = true;
    }

    private Type units = Type.GENERATIONS;

    public final Type getUnits() {
        return units;
    }

    public final void setUnits(Type units) {
        this.units = units;
    }

    private double startTime;
    private double finishTime;
    private int cumulativeSampleCount;
    private int cumulativeCoalescentCount;
    private boolean intervalsKnown;
    private final double[] sampleTimes;
    private final double[] coalescentTimes;
    private final double[] eventTimes;
    private final int eventCount;
    private final int sampleCount;
    private final int coalescentCount;
    private final int intervalCount;

    private final double[] intervals;
    private final int[] lineageCounts;
    private final IntervalType[] intervalTypes;

    /**
     * Testing speed of sorting using Arrays.sort()
     *
     * Take-homes:
     * 1) sorting doubles is nearly 4 times faster than objects
     * 2) parallel sort achieves nothing and may be slower for objects
     * 3) sorting an already sorted array is 15x faster than an unsorted one
     *
     * Times in ms for an array of length 1600:
     * Unsorted doubles, 1000000 reps, time = 12047
     * Unsorted doubles (parallel sort), 1000000 reps, time = 11982
     * Pre-sorted doubles, 1000000 reps, time = 878
     * Unsorted objects, 1000000 reps, time = 45126
     * Unsorted objects (parallel sort), 1000000 reps, time = 50847
     *
     * Test the time taken to do the array copy in each loop:
     * Array copy doubles, 1000000 reps, time = 185
     * Array copy objects, 1000000 reps, time = 80
     * @param args
     */
    public static void main(String[] args) {
        int count = 1600;
        int reps = 1000000;
        double[] randomTimes = new double[count];
        double[] orderedTimes = new double[count];
        double[] times = new double[count];
        Test[] randomTests = new Test[count];
        Test[] tests = new Test[count];

        for (int j = 0; j < count; j++) {
            randomTimes[j] = MathUtils.nextDouble();
            orderedTimes[j] = MathUtils.nextDouble();
            randomTests[j] = new Test(MathUtils.nextDouble());
        }
        Arrays.sort(orderedTimes);

        // Check the time of the array copy in each loop
        long startTime = System.currentTimeMillis();
        for (int i = 0; i < reps; i++) {
            System.arraycopy(randomTimes,0, times, 0, count);
        }
        System.out.println("Array copy doubles, " + reps + " reps, time = " + (System.currentTimeMillis() - startTime));

        startTime = System.currentTimeMillis();
        for (int i = 0; i < reps; i++) {
            System.arraycopy(randomTests,0, tests, 0, count);
        }
        System.out.println("Array copy objects, " + reps + " reps, time = " + (System.currentTimeMillis() - startTime));

        startTime = System.currentTimeMillis();
        for (int i = 0; i < reps; i++) {
            System.arraycopy(randomTimes,0, times, 0, count);
            Arrays.sort(times);
        }
        System.out.println("Unsorted doubles, " + reps + " reps, time = " + (System.currentTimeMillis() - startTime));

        startTime = System.currentTimeMillis();
        for (int i = 0; i < reps; i++) {
            System.arraycopy(randomTimes,0, times, 0, count);
            Arrays.parallelSort(times);
        }
        System.out.println("Unsorted doubles (parallel sort), " + reps + " reps, time = " + (System.currentTimeMillis() - startTime));

        startTime = System.currentTimeMillis();
        for (int i = 0; i < reps; i++) {
            System.arraycopy(orderedTimes,0, times, 0, count);
            Arrays.sort(times);
        }
        System.out.println("Pre-sorted doubles, " + reps + " reps, time = " + (System.currentTimeMillis() - startTime));

        startTime = System.currentTimeMillis();
        for (int i = 0; i < reps; i++) {
            System.arraycopy(orderedTimes,0, times, 0, count);
            for (int j = 0; j < 3; j++) {
                times[MathUtils.nextInt(count)] = MathUtils.nextDouble();
            }
            Arrays.sort(times);
        }
        System.out.println("Mostly sorted doubles, " + reps + " reps, time = " + (System.currentTimeMillis() - startTime));


        startTime = System.currentTimeMillis();
        for (int i = 0; i < reps; i++) {
            System.arraycopy(randomTests,0, tests, 0, count);
            Arrays.sort(tests);
        }
        System.out.println("Unsorted objects, " + reps + " reps, time = " + (System.currentTimeMillis() - startTime));

        startTime = System.currentTimeMillis();
        for (int i = 0; i < reps; i++) {
            System.arraycopy(randomTests,0, tests, 0, count);
            Arrays.parallelSort(tests);
        }
        System.out.println("Unsorted objects (parallel sort), " + reps + " reps, time = " + (System.currentTimeMillis() - startTime));


    }

    static class Test implements Comparable<Test> {
        double value;

        public Test(double value) {
            this.value = value;
        }


        @Override
        public int compareTo(Test o) {
            return Double.compare(value, o.value);
        }
    }
}