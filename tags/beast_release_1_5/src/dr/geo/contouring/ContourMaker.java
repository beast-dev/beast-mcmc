package dr.geo.contouring;

/**
 * @author Marc Suchard
 */

public interface ContourMaker {

    ContourPath[] getContourPaths(double level);

}
