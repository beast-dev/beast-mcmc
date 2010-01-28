package dr.geo;

import dr.math.distributions.MultivariateNormalDistribution;

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
     * @param normal multivariate normal distribution characterizing Brownian bridge
     * @param start starting spacetime
     * @param end   ending spacetime
     * @param depth depth of divide and conquer recursion resulting in 2^(depth-1) new points
     * @param rejector invalid space/time rejector
     * @return list of points, containing 2 + 2^(depth-1) points
     */
    public static List<SpaceTime> divideConquerBrownianBridge(MultivariateNormalDistribution normal, SpaceTime start, SpaceTime end, int depth, SpaceTimeRejector rejector) {

        List<SpaceTime> points = new LinkedList<SpaceTime>();
        points.add(start);
        points.add(end);

        divideConquerBrownianBridge(normal, 0, points, depth, rejector);
        return points;
    }

    public static int divideConquerBrownianBridge(MultivariateNormalDistribution normal, int point0, List<SpaceTime> points, int depth, SpaceTimeRejector rejector) {

        if (depth > 0) {

            SpaceTime pt0 = points.get(point0);
            SpaceTime pt1 = points.get(point0 + 1);

            double t0 = pt0.getTime();
            double[] x0 = pt0.getX();

            double t1 = pt1.getTime();
            double[] x1 = pt1.getX();

            double tm = (t1 + t0) / 2.0;

            double t0m = tm - t1;
            double p0m = 1.0 / t0m;

            double tm1 = t1 - tm;
            double pm1 = 1.0 / tm1;

            // p0m + pm1 = precision of pt at tm,
            double v01 = 1.0 / (p0m + pm1);

            final int dim = x0.length;
            double[] xm = new double[dim];
            for (int i = 0; i < dim; i++) {
                xm[i] = (p0m * x0[i] + pm1 * x1[i]) * v01; 
            }

            SpaceTime s;
            do {
//                getBridgedX(tm, xm, normal, pt0, pt1);
                s = new SpaceTime(tm, normal.nextScaledMultivariateNormal(xm, v01));
            } while (rejector != null && rejector.reject(s));

            points.add(point0 + 1, s);

            int endPoint = divideConquerBrownianBridge(normal, point0, points, depth - 1, rejector);
            return divideConquerBrownianBridge(normal, endPoint, points, depth - 1, rejector);

        } else return point0 + 1;
    }

    //  fill x by bridge between start and end
//    private static void getBridgedX(double t, double[] x, MultivariateNormalDistribution normal, SpaceTime start, SpaceTime end) {
//
//    }

}
