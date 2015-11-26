/*
 * ColourChange.java
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

package dr.evolution.tree;

/**
 * @author Alexei Drummond
 *
 * @version $Id: ColourChange.java,v 1.2 2005/09/21 18:47:58 gerton Exp $
 */
public class ColourChange {

    // time of the change in colour
    private double time;

    // the colour above this time
    private int aboveColour;

    public ColourChange(ColourChange change) {
        this(change.time, change.aboveColour);
    }

    public ColourChange(double time, int aboveColour) {
        this.time = time;
        this.aboveColour = aboveColour;
    }

    /**
     * @return the time of this change in colour
     */
    public final double getTime() { return time; }

    /**
     * @return the colour above the time
     */
    public final int getColourAbove() { return aboveColour; }
    
    /**
     * @param aboveColour sets the colour of this event
     */
    public final void setColourAbove(int aboveColour) { this.aboveColour = aboveColour; }
}

