package dr.geo;

/**
 * @author AlexeiMarc
 */
public class SpaceTime {

    double[] space;
    double time;

    public SpaceTime(double time, double[] space) {
        this.time = time;
        this.space = space;
    }

    public double[] getX() {
        return space;
    }

    public double getTime() {
        return time;
    }
}
