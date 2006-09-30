/*
 * SpaceTime.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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

package dr.app.simcoal;

import dr.evolution.continuous.SphericalPolarCoordinates;

import java.util.Comparator;

/**
 * @author Alexei Drummond
 */
public class SpaceTime {

    double longitude;
    double latitude;
    double time;
    SphericalPolarCoordinates coord;
    String name = null;

    public SpaceTime(double latitude, double longitude, double time) {
        this.longitude = longitude;
        this.latitude = latitude;
        this.time = time;
        coord = new SphericalPolarCoordinates(latitude, longitude);
    }

    public double euclideanDistance(SpaceTime spaceTime) {

        double xDiff = longitude -spaceTime.longitude;
        double yDiff = latitude -spaceTime.latitude;


        return Math.sqrt(xDiff*xDiff + yDiff * yDiff);
    }

    public double distanceInKilometres(SpaceTime spaceTime) {

        double distance = coord.distance(spaceTime.coord);

        return distance;

    }


    public double relativeTime(SpaceTime spaceTime) {
        return time - spaceTime.time;
    }

    public String getName() { return name; }

    public void setName(String name) { this.name = name; }

    public static Comparator<SpaceTime> getLatitudeComparator() {
        return new Comparator<SpaceTime>() {

            public int compare(SpaceTime sp1, SpaceTime sp2) {

                if (sp1.latitude > sp2.latitude) return 1;
                if (sp1.latitude < sp2.latitude) return -1;
                return 0;
            }
        };
    }
}
