package dr.geo;

import java.awt.geom.Point2D;

/**
 * @author Alexei Drummond
 * @author Marc Suchard
 */
public class SpaceTime {

    double[] space;
    double time;

    public SpaceTime(double time, Point2D space) {
        this.time = time;
        this.space = new double[]{space.getX(), space.getY()};
    }

    public SpaceTime(double time, double[] space) {
        this.time = time;
        this.space = space;
    }

    public double[] getX() {
        return space;
    }

    public double getX(int index) {
        return space[index];
    }

    public double getTime() {
        return time;
    }

    public String toString() {
        StringBuilder builder = new StringBuilder();
        builder.append(time);
        for (double s : space) {
            builder.append("\t").append(s);
        }
        return builder.toString();
    }
}
