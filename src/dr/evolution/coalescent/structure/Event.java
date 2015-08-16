/*
 * Event.java
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

/**
 * @author Alexei
 * @version $Id: Event.java,v 1.5 2005/09/21 15:58:34 rambaut Exp $
 */
public class Event implements Comparable {

    double time;
    int[] lineageChanges;
    private IntervalType type;
    private int aboveColour;
    private int belowColour;

    public static Event createCoalescentEvent(double time, int colour, int colourCount) {

        int[] lineageChanges = new int[colourCount];
        lineageChanges[colour] = -1;

        return new Event(time, IntervalType.COALESCENT, lineageChanges, colour, colour);
    }

    public static Event createAddSampleEvent(double time, int colour, int colourCount) {

        int[] lineageChanges = new int[colourCount];
        lineageChanges[colour] = +1;

        return new Event(time, IntervalType.SAMPLE, lineageChanges, colour, colour);
    }

    public static Event createMigrationEvent(double time, int belowColour, int aboveColour, int colourCount) {

        int[] lineageChanges = new int[colourCount];
        lineageChanges[belowColour] = -1;
        lineageChanges[aboveColour] = +1;

        return new Event(time, IntervalType.MIGRATION, lineageChanges, aboveColour, belowColour);
    }

    private Event(double time, IntervalType type, int[] lineageChanges, int aboveColour, int belowColour) {
        this.time = time;
        this.type = type;
        this.lineageChanges = lineageChanges;
        this.aboveColour = aboveColour;
        this.belowColour = belowColour;
    }

    public int getAboveColour() {
        return aboveColour;
    }

    public int getBelowColour() {
        return belowColour;
    }

    public IntervalType getType() {
        return type;
    }

    public String toString() {
        return type + " event at time " + time + " (above=" + aboveColour + ", below=" + belowColour + ")";
    }

    public int compareTo(Object o) {
        Event e2 = (Event) o;

        return new Double(time).compareTo(e2.time);
    }
}
