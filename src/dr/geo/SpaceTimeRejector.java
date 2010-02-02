package dr.geo;

import java.util.List;

/**
 * @author
 */
public interface SpaceTimeRejector {

    boolean reject(SpaceTime point, int attribute);

    // removes all rejects
    void reset();

    List<Reject> getRejects();

//    void setStop(boolean stop);
//
//    boolean getStop();
}
