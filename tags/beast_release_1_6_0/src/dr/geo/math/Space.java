package dr.geo.math;

/**
 * @author Marc Suchard
 */

public enum Space {
    CARTESIAN,  // (x,y,z,...)
    LAT_LONG,   // (lat, long) in S^2
    DEGREES,    // (degrees) in S^1
    RADIANS;    // (radians) in S^1
}
