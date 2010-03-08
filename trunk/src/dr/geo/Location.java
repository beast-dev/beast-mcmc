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
