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
