package dr.geo;

import java.util.LinkedList;
import java.util.List;

/**
 * @author Alexei Drummond
 * @author Marc Suchard
 */
public class MultivariateBrownianBridge {

    /**
     * Divide and conquer brownian bridge
     *
     * @param start starting spacetime
     * @param end   ending spacetime
     * @param depth depth of divide and conquer recursion resulting in 2^(depth-1) new points
     * @return list of points, containing 2 + 2^(depth-1) points
     */
    public static List<SpaceTime> divideConquerBrownianBridge(double[][] precision, SpaceTime start, SpaceTime end, int depth, SpaceTimeRejector rejector) {

        List<SpaceTime> points = new LinkedList<SpaceTime>();
        points.add(start);
        points.add(end);

        divideConquerBrownianBridge(precision, 0, points, depth, rejector);
        return points;
    }

    public static int divideConquerBrownianBridge(double[][] precision, int point0, List<SpaceTime> points, int depth, SpaceTimeRejector rejector) {

        if (depth > 0) {

            SpaceTime p0 = points.get(point0);
            SpaceTime p1 = points.get(point0 + 1);

            double t0 = p0.getTime();
            double[] x0 = p0.getX();

            double t1 = p1.getTime();
            double[] x1 = p1.getX();

            double tm = (t1 - t0) / 2.0 + t0;

            double[] xm = new double[x0.length];
            SpaceTime s;
            do {
                getBridgedX(tm, xm, precision, p0, p1);
                s = new SpaceTime(tm, xm);
            } while (rejector != null && rejector.reject(s));

            points.add(point0 + 1, s);

            int endPoint = divideConquerBrownianBridge(precision, point0, points, depth - 1, rejector);
            return divideConquerBrownianBridge(precision, endPoint, points, depth - 1, rejector);

        } else return point0 + 1;
    }

    //  fill x by bridge between start and end
    private static void getBridgedX(double t, double[] x, double[][] precision, SpaceTime start, SpaceTime end) {

    }

}
