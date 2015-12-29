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

package dr.evolution.util;

import dr.util.Attribute;
import dr.util.Identifiable;

import java.util.TreeMap;
import java.util.Map;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class Location implements Identifiable {

    public static final String LOCATION = "location";

    private Location(final String id, final int index, final String description, final double longitude, final double latitude) {
        this.id = id;
        this.index = index;
        this.description = description;
        this.longitude = longitude;
        this.latitude = latitude;
    }

    public String getId() {
        return id;
    }

    public void setId(String id) {
//        this.id = id;
        // ignore setting...
    }

    public String getDescription() {
        return description;
    }

    public double getLongitude() {
        return longitude;
    }

    public double getLatitude() {
        return latitude;
    }

    public int getIndex() {
        return index;
    }

    @Override
    public String toString() {
        return getId();
    }

    private final String id;
    private final int index;
    private final String description;

    private final double longitude;
    private final double latitude;

    // STATIC MEMBERS

    public static Location newLocation(String id, String description, double longitude, double latitude) {
        if (getLocation(id) != null) {
            throw new IllegalArgumentException("Location with id, " + id + ", already exists");
        }
        int index = getLocationCount();
        Location location =  new Location(id, index, description, longitude, latitude);
        locations.put(index, location);
        return location;
    }

    public static Location getLocation(String id) {
        for (Location location : locations.values()) {
            if (location.getId().equals(id)) {
                return location;
            }
        }
        return null;
    }

    public static int getLocationCount() {
        return locations.keySet().size();
    }

    public static Location getLocation(int index) {
        return locations.get(index);
    }

    static private Map<Integer, Location> locations = new TreeMap<Integer, Location>();
}
