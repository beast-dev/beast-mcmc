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

import dr.evolution.util.Units;

import java.util.Arrays;

/**
 * A concrete class for a set of coalescent intervals.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 */
public class Intervals implements IntervalList {


    public Intervals(int maxEventCount,boolean eventsNeedSorting) {
        startTime = Double.POSITIVE_INFINITY;

        events = new Event[maxEventCount];
        for (int i = 0; i < maxEventCount; i++) {
            events[i] = new Event();
        }
        eventCount = 0;
        sampleCount = 0;

        intervals = new double[maxEventCount - 1];
        intervalTypes = new IntervalType[maxEventCount - 1];
        lineageCounts = new int[maxEventCount - 1];

        intervalsKnown = false;
        this.eventsNeedSorting=eventsNeedSorting;
    }
    public Intervals(int maxEventCount) {
        this(maxEventCount, true);
    }

    public void copyIntervals(Intervals source) {
        intervalsKnown = source.intervalsKnown;
        eventCount = source.eventCount;
        sampleCount = source.sampleCount;
        intervalCount = source.intervalCount;
        startTime = source.startTime;

        //don't copy the actual events..
        /*
          for (int i = 0; i < events.length; i++) {
              events[i].time = source.events[i].time;
              events[i].type = source.events[i].type;
          }*/

        if (intervalsKnown) {
            System.arraycopy(source.intervals, 0, intervals, 0, intervals.length);
            System.arraycopy(source.intervalTypes, 0, intervalTypes, 0, intervals.length);
            System.arraycopy(source.lineageCounts, 0, lineageCounts, 0, intervals.length);
        }
    }

    public void resetEvents() {
        startTime = Double.POSITIVE_INFINITY;

        intervalsKnown = false;
        eventCount = 0;
        sampleCount = 0;
    }

    public void addSampleEvent(double time) {
        this.addSampleEvent(time,-1);
    }
    public void addSampleEvent(double time,int nodeNumber) {
        if (time < startTime) {
            startTime = time;
        }

        events[eventCount].time = time;
        events[eventCount].type = IntervalType.SAMPLE;
        events[eventCount].nodeNumber=nodeNumber;
        eventCount++;
        sampleCount++;
        intervalsKnown = false;
    }
    public void addCoalescentEvent(double time) {
        this.addCoalescentEvent(time, -1);
    }

    public void addCoalescentEvent(double time, int nodeNumber) {
        events[eventCount].time = time;
        events[eventCount].type = IntervalType.COALESCENT;
        events[eventCount].nodeNumber=nodeNumber;
        eventCount++;
        intervalsKnown = false;
    }

    public void addMigrationEvent(double time, int destination) {
        events[eventCount].time = time;
        events[eventCount].type = IntervalType.MIGRATION;
        events[eventCount].info = destination;
        events[eventCount].nodeNumber=-1;
        eventCount++;
        intervalsKnown = false;
    }

    public void addNothingEvent(double time) {
        events[eventCount].time = time;
        events[eventCount].type = IntervalType.NOTHING;
        events[eventCount].nodeNumber = -1;
        eventCount++;
        intervalsKnown = false;
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
        return events[i].time;
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
    //Return the node that triggers the event
    public int getNodeForEvent(int i){
        if (!intervalsKnown){
            calculateIntervals();
        }
        return events[i].nodeNumber;
    }
    public double getTotalDuration() {

        if (!intervalsKnown) {
            calculateIntervals();
        }
        return events[eventCount - 1].time;
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

        if(eventsNeedSorting) {
            Arrays.sort(events, 0, eventCount);
        }
        if (events[0].type != IntervalType.SAMPLE) {
            throw new IllegalArgumentException("First event is not a sample event");
        }

        intervalCount = eventCount - 1;

        double lastTime = events[0].time;

        int lineages = 1;
        for (int i = 1; i < eventCount; i++) {

            intervals[i - 1] = events[i].time - lastTime;
            intervalTypes[i - 1] = events[i].type;
            lineageCounts[i - 1] = lineages;
            if (events[i].type == IntervalType.SAMPLE) {
                lineages++;
            } else if (events[i].type == IntervalType.COALESCENT) {
                lineages--;
            }
            lastTime = events[i].time;
        }
        intervalsKnown = true;
    }

    private Units.Type units = Units.Type.GENERATIONS;

    public final Units.Type getUnits() {
        return units;
    }

    public final void setUnits(Type units) {
        this.units = units;
    }

    private class Event implements Comparable {

        public int compareTo(Object o) {
            double t = ((Event) o).time;
            if (t < time) {
                return 1;
            } else if (t > time) {
                return -1;
            } else {
                // events are at exact same time so sort by type
                return type.compareTo(((Event) o).type);
            }
        }

        /**
         * The type of event
         */
        IntervalType type;

        /**
         * The time of the event
         */
        double time;

        /**
         * Some extra information for the event (e.g., destination of a migration)
         */
        int info;
        int nodeNumber;

    }

    private double startTime;

    private Event[] events;
    private int eventCount;
    private int sampleCount;

    private boolean intervalsKnown = false;
    private final boolean eventsNeedSorting;
    private double[] intervals;
    private int[] lineageCounts;
    private IntervalType[] intervalTypes;
    private int intervalCount = 0;
}