package dr.geo;

/**
 * Created by IntelliJ IDEA.
 * User: adru001
 * Date: Feb 2, 2010
 * Time: 11:27:54 AM
 * To change this template use File | Settings | File Templates.
 */
public class Reject {

    private int depth;
    private SpaceTime spaceTime;

    public Reject(int depth, SpaceTime s) {
        this.depth = depth;
        this.spaceTime = s;
    }

    public int getDepth() {
        return depth;
    }

    public SpaceTime getSpaceTime() {
        return spaceTime;
    }
}
