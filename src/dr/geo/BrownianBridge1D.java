/*
 * BrownianBridge1D.java
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

import dr.math.MathUtils;

import java.awt.geom.Point2D;
import java.util.LinkedList;
import java.util.List;

/**
 * @author Alexei Drummond
 */

public class BrownianBridge1D {

    /**
     * Divide and conquer brownian bridge
     *
     * @param t0    starting time
     * @param t1    end time
     * @param y0    starting value
     * @param y1    ending value
     * @param depth depth of divide and conquer recursion resulting in 2^(depth-1) new points
     * @return list of points, containing 2 + 2^(depth-1) points
     */
    public static List<Point2D> divideConquerBrownianBridge(double D, double t0, double y0, double t1, double y1, int depth, SpaceTimeRejector1D rejector) {

        return divideConquerBrownianBridge(D, new Point2D.Double(t0, y0), new Point2D.Double(t1, y1), depth, rejector);
    }


    /**
     * Divide and conquer brownian bridge
     *
     * @param t0y0  starting point in time and 1d space
     * @param t1y1  end point in time and 1d space
     * @param depth depth of divide and conquer recursion resulting in 2^(depth-1) new points
     * @return list of points, containing 2 + 2^(depth-1) points
     */
    public static List<Point2D> divideConquerBrownianBridge(double D, Point2D t0y0, Point2D t1y1, int depth, SpaceTimeRejector1D rejector) {

        List<Point2D> points = new LinkedList<Point2D>();
        points.add(t0y0);
        points.add(t1y1);

        divideConquerBrownianBridge(D, 0, points, depth, rejector);
        return points;
    }

    public static int divideConquerBrownianBridge(double D, int point0, List<Point2D> points, int depth, SpaceTimeRejector1D rejector) {

        if (depth > 0) {

            Point2D p0 = points.get(point0);
            Point2D p1 = points.get(point0 + 1);

            double t0 = p0.getX();
            double y0 = p0.getY();

            double t1 = p1.getX();
            double y1 = p1.getY();

            double tm = (t1 - t0) / 2.0 + t0;

            double mean = y0 + 0.5 * (y1 - y0);
            double stdev = Math.sqrt(D * (t1 - t0) / 4.0);

            double xm = MathUtils.nextGaussian() * stdev + mean;
            if (rejector != null) {
                while (rejector.reject(tm, xm)) {
                    xm = MathUtils.nextGaussian() * stdev + mean;
                }
            }

            points.add(point0 + 1, new Point2D.Double(tm, xm));

            int endPoint = divideConquerBrownianBridge(D, point0, points, depth - 1, rejector);
            return divideConquerBrownianBridge(D, endPoint, points, depth - 1, rejector);

        } else return point0 + 1;
    }

    public static void main(String[] args) {

        List<Point2D> points = divideConquerBrownianBridge(1, 0, 0, 4, 0, 5, null);

        for (Point2D p : points) {
            System.out.println(p.getX() + "\t" + p.getY());
        }
    }

}
