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

import java.util.Arrays;

/**
 * A concrete class for a set of coalescent intervals.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
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
        throw new UnsupportedOperationException("not supported in FastIntervals");
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
        throw new UnsupportedOperationException("not supported in FastIntervals");
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
        return startTime;
    }

    public IntervalType getIntervalType(int i) {
        if (!intervalsKnown) {
            calculateIntervals();
        }
        return intervalTypes[i];
    }
    public double getTotalDuration() {
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

//        System.arraycopy(sampleTimes, 0, eventTimes, 0, sampleCount);
//        System.arraycopy(coalescentTimes, 0, eventTimes, sampleCount, coalescentCount);
        Arrays.sort(sampleTimes);
        Arrays.sort(coalescentTimes);

        // sample should come first
        double lastTime = sampleTimes[0];

        int s = 1;
        int c = 0;
        int i = 0;
        int lineages = 1;

        while (s < sampleCount && c < coalescentCount) {
            if (sampleTimes[s] <= coalescentTimes[c]) {
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
    private int cumulativeSampleCount = 0;
    private int cumulativeCoalescentCount = 0;
    private boolean intervalsKnown = false;

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
}