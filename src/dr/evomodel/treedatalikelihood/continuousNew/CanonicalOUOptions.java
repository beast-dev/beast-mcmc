package dr.evomodel.treedatalikelihood.continuousNew;

import java.util.Locale;

/**
 * Parser-facing options for the canonical OU tree backend.
 */
public final class CanonicalOUOptions {

    public static final String COORDINATES_ATTRIBUTE = "canonicalOUCoordinates";
    public static final String BLOCK_BASIS_HOMOGENEOUS_PROPERTY =
            "beast.experimental.canonicalOU.blockBasisHomogeneous";

    public static final String COORDINATES_AUTO = "auto";
    public static final String COORDINATES_NATURAL = "natural";
    public static final String COORDINATES_HOMOGENEOUS_BLOCK_BASIS = "homogeneousBlockBasis";

    public static final String[] ALLOWED_COORDINATES = {
            COORDINATES_AUTO,
            COORDINATES_NATURAL,
            COORDINATES_HOMOGENEOUS_BLOCK_BASIS
    };

    public enum Coordinates {
        AUTO,
        NATURAL,
        HOMOGENEOUS_BLOCK_BASIS
    }

    private static final CanonicalOUOptions AUTO_OPTIONS =
            new CanonicalOUOptions(Coordinates.AUTO);

    private final Coordinates coordinates;

    private CanonicalOUOptions(final Coordinates coordinates) {
        this.coordinates = coordinates;
    }

    public static CanonicalOUOptions defaults() {
        return AUTO_OPTIONS;
    }

    public static CanonicalOUOptions natural() {
        return new CanonicalOUOptions(Coordinates.NATURAL);
    }

    public static CanonicalOUOptions homogeneousBlockBasis() {
        return new CanonicalOUOptions(Coordinates.HOMOGENEOUS_BLOCK_BASIS);
    }

    public static CanonicalOUOptions parseCoordinates(final String rawCoordinates) {
        final String coordinates = rawCoordinates.trim().toLowerCase(Locale.ROOT);
        if (COORDINATES_AUTO.toLowerCase(Locale.ROOT).equals(coordinates)) {
            return defaults();
        }
        if (COORDINATES_NATURAL.toLowerCase(Locale.ROOT).equals(coordinates)) {
            return natural();
        }
        if (COORDINATES_HOMOGENEOUS_BLOCK_BASIS.toLowerCase(Locale.ROOT).equals(coordinates)) {
            return homogeneousBlockBasis();
        }
        throw new IllegalArgumentException(
                COORDINATES_ATTRIBUTE + " must be one of "
                        + "\"" + COORDINATES_AUTO + "\", "
                        + "\"" + COORDINATES_NATURAL + "\", or "
                        + "\"" + COORDINATES_HOMOGENEOUS_BLOCK_BASIS + "\", not \""
                        + rawCoordinates + "\".");
    }

    public boolean isAutomatic() {
        return coordinates == Coordinates.AUTO;
    }

    public boolean isExplicit() {
        return coordinates != Coordinates.AUTO;
    }

    public boolean useHomogeneousBlockBasis() {
        switch (coordinates) {
            case HOMOGENEOUS_BLOCK_BASIS:
                return true;
            case NATURAL:
                return false;
            case AUTO:
                return Boolean.parseBoolean(
                        System.getProperty(BLOCK_BASIS_HOMOGENEOUS_PROPERTY, "false"));
            default:
                throw new IllegalStateException("Unhandled canonical OU coordinates: " + coordinates);
        }
    }

    public String describeSelectionSource() {
        if (coordinates == Coordinates.AUTO) {
            return BLOCK_BASIS_HOMOGENEOUS_PROPERTY;
        }
        return COORDINATES_ATTRIBUTE + "=\"" + toXmlString() + "\"";
    }

    public String toXmlString() {
        switch (coordinates) {
            case AUTO:
                return COORDINATES_AUTO;
            case NATURAL:
                return COORDINATES_NATURAL;
            case HOMOGENEOUS_BLOCK_BASIS:
                return COORDINATES_HOMOGENEOUS_BLOCK_BASIS;
            default:
                throw new IllegalStateException("Unhandled canonical OU coordinates: " + coordinates);
        }
    }
}
