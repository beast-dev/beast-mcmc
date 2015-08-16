/*
 * SpaceTime.java
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

import java.awt.*;
import java.awt.geom.AffineTransform;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Point2D;

/**
 * @author Alexei Drummond
 * @author Marc Suchard
 */
public class SpaceTime {

    double[] space;
    double time;

    public SpaceTime(SpaceTime s) {

        this.space = new double[s.space.length];
        System.arraycopy(s.space, 0, this.space, 0, s.space.length);
        this.time = s.time;
    }

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

    public static void paintDot(SpaceTime s, double radius, AffineTransform transform, Graphics2D g2d) {

        Point2D pointRaw = new Point2D.Double(s.getX(0), s.getX(1));
        Point2D pointT = new Point2D.Double();

        transform.transform(pointRaw, pointT);

        Shape pointShape = new Ellipse2D.Double(pointT.getX() - radius, pointT.getY() - radius, 2.0 * radius, 2.0 * radius);

        g2d.fill(pointShape);
    }

}
