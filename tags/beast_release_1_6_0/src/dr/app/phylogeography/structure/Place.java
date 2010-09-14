package dr.app.phylogeography.structure;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class Place extends GeoItem {
    public Place(final String name, final String description,
                 final Coordinates coordinates,
                 final double startTime, final double duration) {
        super(name, startTime, duration);
        this.description = description;
        this.coordinates = coordinates;
    }

    public String getDescription() {
        return description;
    }

    public Coordinates getCoordinates() {
        return coordinates;
    }

    private final String description;
    private final Coordinates coordinates;
}
