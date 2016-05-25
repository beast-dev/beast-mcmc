/*
 * Location.java
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

package dr.geo;

/**
 * @author Alexei Drummond
 */
public class Location {

    Location(int i, int j) {
        this.i = i;
        this.j = j;
    }

    public final int i, j;

    public boolean equals(Object obj) {
        if (!(obj instanceof Location)) return false;
        Location loc = (Location) obj;

        return loc.i == i && loc.j == j;
    }

    public String toString() {
        return "(" + i + ", " + j + ")";
    }

    public int hashCode() {
        return (i * 1009 + j);
    }

}
