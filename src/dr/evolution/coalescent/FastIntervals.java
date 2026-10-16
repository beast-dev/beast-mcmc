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
        sortedSampleTimes = new double[sampleCount];
        sortedCoalescentTimes = new double[coalescentCount];
        sampleSorter = new IncrementalSorter(sampleCount);
        coalescentSorter = new IncrementalSorter(coalescentCount);
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

        // keep the sort orders current, so the copy can update incrementally if it is swapped back in
        sampleSorter.copy(source.sampleSorter);
        coalescentSorter.copy(source.coalescentSorter);
    }

    @Override
    public void resetEvents() {
        startTime = Double.POSITIVE_INFINITY;
        finishTime = Double.NEGATIVE_INFINITY;

        intervalsKnown = false;
        cumulativeSampleCount = 0;
        cumulativeCoalescentCount = 0;
        eventsInNodeOrder = true;
    }

    @Override
    public void addSampleEvent(double time) {
        // without node numbers the events may come in any order
        eventsInNodeOrder = false;
        addSample(time);
    }

    private void addSample(double time) {
        if (time < startTime) {
            startTime = time;
        }

        sampleTimes[cumulativeSampleCount] = time;
        cumulativeSampleCount++;
        intervalsKnown = false;
    }

    @Override
    public void addSampleEvent(double time, int nodeNumber) {
        // events added with node numbers come in the same order each time (see calculateIntervals)
        addSample(time);
    }

    @Override
    public void addCoalescentEvent(double time) {
        // without node numbers the events may come in any order
        eventsInNodeOrder = false;
        addCoalescent(time);
    }

    private void addCoalescent(double time) {
        if (time > finishTime) {
            finishTime = time;
        }

        coalescentTimes[cumulativeCoalescentCount] = time;
        cumulativeCoalescentCount++;
        intervalsKnown = false;
    }

    @Override
    public void addCoalescentEvent(double time, int nodeNumber) {
        // events added with node numbers come in the same order each time (see calculateIntervals)
        addCoalescent(time);
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

        if (eventsInNodeOrder) {
            // Events added with node numbers arrive in the same order every time (TreeIntervals adds them
            // node by node), so only the times that have changed need to be placed in the previous order.
            sampleSorter.sort(sampleTimes, sortedSampleTimes);
            coalescentSorter.sort(coalescentTimes, sortedCoalescentTimes);
        } else {
            System.arraycopy(sampleTimes, 0, sortedSampleTimes, 0, sampleCount);
            System.arraycopy(coalescentTimes, 0, sortedCoalescentTimes, 0, coalescentCount);
            Arrays.sort(sortedSampleTimes);
            Arrays.sort(sortedCoalescentTimes);
            sampleSorter.invalidate();
            coalescentSorter.invalidate();
        }
        final double[] sampleTimes = sortedSampleTimes;
        final double[] coalescentTimes = sortedCoalescentTimes;

        // sample should come first
        double lastTime = sampleTimes[0];

        int s = 1;
        int c = 0;
        int i = 0;
        int lineages = 1;

//         eventTimes[0] = sampleTimes[0]; // keeping in case we need to revert

        while (s < sampleCount || c < coalescentCount) {
            if (s < sampleCount && sampleTimes[s] <= coalescentTimes[c]) {
                intervals[i] = sampleTimes[s] - lastTime;
                eventTimes[i] = sampleTimes[s];                
//                 eventTimes[i + 1] = sampleTimes[s];
                intervalTypes[i] = IntervalType.SAMPLE;
                lineageCounts[i] = lineages;
                lastTime = sampleTimes[s];
                lineages++;
                s++;
            } else {
                intervals[i] = coalescentTimes[c] - lastTime;
				eventTimes[i] = coalescentTimes[c];                
//                 eventTimes[i + 1] = coalescentTimes[c];
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
    private final double[] sortedSampleTimes;
    private final double[] sortedCoalescentTimes;
    private final IncrementalSorter sampleSorter;
    private final IncrementalSorter coalescentSorter;
    private boolean eventsInNodeOrder = true;
    private final double[] eventTimes;
    private final int eventCount;
    private final int sampleCount;
    private final int coalescentCount;
    private final int intervalCount;

    private final double[] intervals;
    private final int[] lineageCounts;
    private final IntervalType[] intervalTypes;

    /**
     * Sorts an array of times that arrive in the same order every time (one per node), where usually
     * only a few have changed since the last sort. It keeps the order (as indices into the array) found
     * by the last sort: the unchanged times keep their relative order, so only the changed ones need
     * sorting, after which the two are merged. This takes O(n + k log k) for k changed times, instead
     * of O(n log n), and gives exactly the same sorted values as Arrays.sort.
     */
    static final class IncrementalSorter {
        private final double[] previous;
        private int[] order;
        private int[] newOrder;
        private final int[] changed;
        private final boolean[] isChanged;
        private final int[] scratch;
        private boolean valid = false;

        IncrementalSorter(int n) {
            previous = new double[n];
            order = new int[n];
            newOrder = new int[n];
            changed = new int[n];
            isChanged = new boolean[n];
            scratch = new int[n];
        }

        void invalidate() {
            valid = false;
        }

        void copy(IncrementalSorter source) {
            valid = source.valid;
            if (valid) {
                System.arraycopy(source.previous, 0, previous, 0, previous.length);
                System.arraycopy(source.order, 0, order, 0, order.length);
            }
        }

        /**
         * @param times  the times in their arrival order
         * @param sorted filled with the times in increasing order
         */
        void sort(double[] times, double[] sorted) {
            final int n = times.length;
            if (!valid) {
                for (int i = 0; i < n; i++) {
                    order[i] = i;
                }
                sortIndices(order, n, times, scratch);
                valid = true;
            } else {
                int k = 0;
                for (int i = 0; i < n; i++) {
                    if (times[i] != previous[i]) {
                        changed[k++] = i;
                        isChanged[i] = true;
                    }
                }
                if (k > 0) {
                    sortIndices(changed, k, times, scratch);
                    // merge the unchanged indices, in their previous order, with the sorted changed ones
                    int r = 0;
                    int c = 0;
                    for (int j = 0; j < n; j++) {
                        int i = order[j];
                        if (isChanged[i]) {
                            continue;
                        }
                        while (c < k && times[changed[c]] < times[i]) {
                            newOrder[r++] = changed[c++];
                        }
                        newOrder[r++] = i;
                    }
                    while (c < k) {
                        newOrder[r++] = changed[c++];
                    }
                    for (int j = 0; j < k; j++) {
                        isChanged[changed[j]] = false;
                    }
                    int[] tmp = order;
                    order = newOrder;
                    newOrder = tmp;
                }
            }
            System.arraycopy(times, 0, previous, 0, n);
            for (int j = 0; j < n; j++) {
                sorted[j] = times[order[j]];
            }
        }

        /**
         * Stable merge sort of the first n indices by their times.
         */
        private static void sortIndices(int[] indices, int n, double[] times, int[] scratch) {
            if (n < 2) {
                return;
            }
            if (n <= 16) {
                for (int i = 1; i < n; i++) {
                    int index = indices[i];
                    double time = times[index];
                    int j = i - 1;
                    while (j >= 0 && times[indices[j]] > time) {
                        indices[j + 1] = indices[j];
                        j--;
                    }
                    indices[j + 1] = index;
                }
                return;
            }
            int[] from = indices;
            int[] to = scratch;
            for (int width = 1; width < n; width *= 2) {
                for (int left = 0; left < n; left += 2 * width) {
                    int middle = Math.min(left + width, n);
                    int right = Math.min(left + 2 * width, n);
                    int i = left;
                    int j = middle;
                    int r = left;
                    while (i < middle && j < right) {
                        to[r++] = times[from[j]] < times[from[i]] ? from[j++] : from[i++];
                    }
                    while (i < middle) {
                        to[r++] = from[i++];
                    }
                    while (j < right) {
                        to[r++] = from[j++];
                    }
                }
                int[] tmp = from;
                from = to;
                to = tmp;
            }
            if (from != indices) {
                System.arraycopy(from, 0, indices, 0, n);
            }
        }
    }

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