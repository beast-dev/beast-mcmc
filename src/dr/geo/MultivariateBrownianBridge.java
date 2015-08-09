/*
 * MultivariateBrownianBridge.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.geo;

import dr.math.distributions.MultivariateNormalDistribution;

import java.util.LinkedList;
import java.util.List;

/**
 * @author Alexei Drummond
 * @author Marc Suchard
 */
public class MultivariateBrownianBridge {

    private static boolean stop;

    /**
     * Divide and conquer brownian bridge; NOT THREAD-SAFE!!!!!
     *
     * @param normal   multivariate normal distribution characterizing Brownian bridge
     * @param start    starting spacetime
     * @param end      ending spacetime
     * @param depth    depth of divide and conquer recursion resulting in 2^(depth-1) new points
     * @param maxTries maximum number of rejections for each level - MUST BE 1 for correct sampling under absorptive boundaries (says AlexeiMarc)
     * @param rejector invalid space/time rejector
     * @return list of points, containing 2 + 2^(depth-1) points
     */
    public static List<SpaceTime> divideConquerBrownianBridge(MultivariateNormalDistribution normal,
                                                              SpaceTime start, SpaceTime end,
                                                              int depth, int maxTries, SpaceTimeRejector rejector) {

        List<SpaceTime> points = new LinkedList<SpaceTime>();
        points.add(start);
        points.add(end);

        stop = false;
        if (divideConquerBrownianBridge(normal, 0, points, depth, maxTries, rejector) == (2 << (depth - 1))) {
            return points;
        }
        stop = false;
        return null;
    }

    public static int divideConquerBrownianBridge(MultivariateNormalDistribution normal,
                                                  int point0, List<SpaceTime> points,
                                                  int depth, int maxTries, SpaceTimeRejector rejector) {

        if (depth > 0 && !stop) {

            SpaceTime pt0 = points.get(point0);
            SpaceTime pt1 = points.get(point0 + 1);

            double t0 = pt0.getTime();
            double[] x0 = pt0.getX();

            double t1 = pt1.getTime();
            double[] x1 = pt1.getX();

            double tm = (t1 + t0) / 2.0;

            double t0m = tm - t0;
            double p0m = 1.0 / t0m;

            double tm1 = t1 - tm;
            double pm1 = 1.0 / tm1;

            // p0m + pm1 = precision of pt at tm,
            double v01 = 1.0 / (p0m + pm1);


//            System.err.println("p0m = " + p0m);
//            System.err.println("pm1 = " + pm1);
//            System.err.println("v01 = " + v01);

            final int dim = x0.length;
            double[] xm = new double[dim];
            for (int i = 0; i < dim; i++) {
                xm[i] = (p0m * x0[i] + pm1 * x1[i]) * v01;
            }

//            System.err.println("Mean: " + new Vector(xm));


            int tries = 0;

            SpaceTime s;
            do {
                s = new SpaceTime(tm, normal.nextScaledMultivariateNormal(xm, v01));
                tries += 1;
                if (tries > maxTries) {
                    stop = true;
                    return 0;
                }
            } while (rejector != null && rejector.reject(s.time, s.space));

            points.add(point0 + 1, s);

            int endPoint = divideConquerBrownianBridge(normal, point0, points, depth - 1, maxTries, rejector);
            return divideConquerBrownianBridge(normal, endPoint, points, depth - 1, maxTries, rejector);

        } else return point0 + 1;
    }
}
