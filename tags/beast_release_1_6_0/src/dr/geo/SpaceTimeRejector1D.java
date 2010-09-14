package dr.geo;

/**
 * @author Alexei Drummond
 */
public interface SpaceTimeRejector1D {

    boolean reject(double time, double y);
}
