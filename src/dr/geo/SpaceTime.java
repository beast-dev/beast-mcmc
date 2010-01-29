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
