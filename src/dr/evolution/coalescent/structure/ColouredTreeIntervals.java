/*
 * ColouredTreeIntervals.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
 */

package dr.evolution.coalescent.structure;

import dr.evolution.coalescent.IntervalType;
import dr.evolution.colouring.BranchColouring;
import dr.evolution.colouring.TreeColouring;
import dr.evolution.tree.ColourChange;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * @author Alexei Drummond
 * @version $Id: ColouredTreeIntervals.java,v 1.9 2005/12/08 13:48:41 rambaut Exp $
 */
public class ColouredTreeIntervals implements StructuredIntervalList {

    int intervalCount;
    int sampleCount;
    List<Event> eventList = new ArrayList<Event>();
    int colourStateCount;
    int[][] lineageCount;
    Type units;

    public ColouredTreeIntervals(Tree tree, TreeColouring colouring) {

        colourStateCount = colouring.getColourCount();
        units = tree.getUnits();
        extractCoalescentEvents(tree, colouring, eventList);
        sampleCount = extractSampleEvents(tree, colouring, eventList);
        extractMigrationEvents(tree, colouring, eventList);
        Collections.sort(eventList);

        lineageCount = new int[eventList.size() + 1][colourStateCount];
        int externalNodeCount = tree.getExternalNodeCount();
        for (int i = 0; i < externalNodeCount; i++) {
            NodeRef node = tree.getExternalNode(i);
            double time = tree.getNodeHeight(node);
            int colour = colouring.getNodeColour(node);
            if (time == 0.0) {
                lineageCount[0][colour] += 1;
            }
        }
        for (int i = 0; i < eventList.size(); i++) {
            Event event = eventList.get(i);
            //System.out.println(event);

            for (int j = 0; j < colourStateCount; j++) {
                lineageCount[i + 1][j] = lineageCount[i][j] + event.lineageChanges[j];
                if (lineageCount[i + 1][j] < 0) {
                    throw new RuntimeException("lineageCount[" + (i + 1) + "][" + j + "] = " + lineageCount[i + 1][j] + ". This is wrong!");
                }
            }
        }
        intervalCount = eventList.size();
    }

    public int getPopulationCount() {
        return lineageCount[0].length;
    }

    /**
     * get number of intervals
     */
    public int getIntervalCount() {
        return intervalCount;
    }

    /**
     * get the total number of sampling events.
     */
    public int getSampleCount() {
        return sampleCount;
    }

    /**
     * Gets an interval.
     */
    public double getInterval(int i) {
        if (i == 0) return getEvent(i).time;
        return getEvent(i).time - getEvent(i - 1).time;
    }

    /**
     * @param population the population of interest
     * @return the number of lineages residing in the given population over this interval.
     */
    public int getLineageCount(int interval, int population) {
        return lineageCount[interval][population];
    }

    /**
     * @param interval the interval of interest
     * @return the number of lineages residing in all populations over this interval.
     */
    public int getLineageCount(int interval) {
        int totalLineages = 0;
        for (int i = 0; i < colourStateCount; i++) {
            totalLineages += lineageCount[interval][i];
        }
        return totalLineages;
    }

    /**
     * Returns the number coalescent events in an interval
     */
    public int getCoalescentEvents(int i) {
        if (getEvent(i).getType() == IntervalType.COALESCENT) return 1;
        return 0;
    }

    /**
     * Returns the type of interval observed.
     */
    public IntervalType getIntervalType(int i) {
        return getEvent(i).getType();
    }

    /**
     * get the total duration of these intervals.
     */
    public double getTotalDuration() {
        return getEvent(getIntervalCount()).time;
    }

    /**
     * Checks whether this set of coalescent intervals is fully resolved
     * (i.e. whether is has exactly one coalescent event in each
     * subsequent interval)
     */
    public boolean isBinaryCoalescent() {
        return true;
    }

    /**
     * Checks whether this set of coalescent intervals coalescent only
     * (i.e. whether is has exactly one or more coalescent event in each
     * subsequent interval)
     */
    public boolean isCoalescentOnly() {
        return false;
    }

    /**
     * Gets the units for this object.
     */
    public Type getUnits() {
        return units;
    }

    /**
     * Sets the units for this object.
     */
    public void setUnits(Type units) {
        this.units = units;
    }

    private void extractCoalescentEvents(Tree tree, TreeColouring colouring, List<Event> eventList) {
        int internalNodeCount = tree.getInternalNodeCount();
        for (int i = 0; i < internalNodeCount; i++) {
            NodeRef node = tree.getInternalNode(i);
            double time = tree.getNodeHeight(node);
            int colour = colouring.getNodeColour(node);
            eventList.add(Event.createCoalescentEvent(time, colour, colouring.getColourCount()));
        }
    }

    private int extractSampleEvents(Tree tree, TreeColouring colouring, List<Event> eventList) {
        int externalNodeCount = tree.getExternalNodeCount();
        int sampleEventCount = 0;
        for (int i = 0; i < externalNodeCount; i++) {
            NodeRef node = tree.getExternalNode(i);
            double time = tree.getNodeHeight(node);
            int colour = colouring.getNodeColour(node);
            if (time != 0.0) {
                eventList.add(Event.createAddSampleEvent(time, colour, colouring.getColourCount()));
                sampleEventCount += 1;
            }
        }
        return sampleEventCount;
    }

    private void extractMigrationEvents(Tree tree, TreeColouring colouring, List<Event> eventList) {
        int nodeCount = tree.getNodeCount();
        for (int i = 0; i < nodeCount; i++) {
            NodeRef node = tree.getNode(i);
            if (!tree.isRoot(node)) {

                BranchColouring branchColouring = colouring.getBranchColouring(node);
                List<ColourChange> changes = branchColouring.getColourChanges();
                int belowColour = colouring.getNodeColour(node);
                for (ColourChange change : changes) {

                    double time = change.getTime();
                    int aboveColour = change.getColourAbove();
                    eventList.add(Event.createMigrationEvent(time, belowColour, aboveColour, colouring.getColourCount()));
                    belowColour = aboveColour;
                }
            }
        }
    }

    public final Event getEvent(int i) {
        return eventList.get(i);
    }

}
