/*
 * Polygon2D.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.xml.XMLParser;
import org.jdom.Element;

import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.util.ArrayList;
import java.util.List;


/**
 * @author Marc A. Suchard
 * @author Philippe Lemey
 */

public class Polygon2D extends AbstractPolygon2D {

    public Polygon2D(double[] x, double[] y) {
        if (x.length != y.length) {
            throw new RuntimeException("Unbalanced arrays");
        }

        if (x[0] != x[x.length - 1] && y[0] != y[y.length - 1]) {
            double[] newX = new double[x.length + 1];
            double[] newY = new double[y.length + 1];
            System.arraycopy(x, 0, newX, 0, x.length);
            System.arraycopy(y, 0, newY, 0, y.length);
            newX[x.length] = x[0];
            newY[y.length] = y[0];
            this.x = newX;
            this.y = newY;
        } else {
            this.x = x;
            this.y = y;
        }
        length = this.x.length - 1;

    }

    public Polygon2D(List<Point2D> points, boolean closed) {
        this.point2Ds = points;
        if (!closed) {
            Point2D start = points.get(0);
            points.add(start);
        }
        convertPointsToArrays();
        length = points.size() - 1;
    }

    public Polygon2D() {
        length = 0;
        point2Ds = new ArrayList<Point2D>();
    }

    public Polygon2D(Element e) {
        List<Element> children = e.getChildren();
        id = e.getAttributeValue(XMLParser.ID);

        parseCoordinates(e);
    }

    public double getProbability(Point2D Point2D, boolean outside) {
        boolean contains = containsPoint2D(Point2D);
        if (outside ^ contains) {
            return 1.0;
        } else {
            return 0.0;
        }
    }

    public double getLogProbability(Point2D Point2D, boolean outside) {
        boolean contains = containsPoint2D(Point2D);
        if (outside ^ contains) {
            return 0.0;
        } else {
            return Double.NEGATIVE_INFINITY;
        }
    }

    @Override
    public void setFillValue(double value) {
        throw new RuntimeException("setFillValue() call not allowed in Polygon2D; try Polygon2DSampling.");
    }

    @Override
    public double getFillValue() {
        throw new RuntimeException("getFillValue() call not allowed in Polygon2D; try Polygon2DSampling.");
    }

    @Override
    public double getLogFillValue() {
        throw new RuntimeException("getLogFillValue() call not allowed in Polygon2D; try Polygon2DSampling.");
    }

    @Override
    public boolean hasFillValue() {
        return false;
    }

    public static void main(String[] args) {

        Polygon2D polygon = new Polygon2D();
        polygon.addPoint2D(new Point2D.Double(-10, -10));
        polygon.addPoint2D(new Point2D.Double(-10, 50));
        polygon.addPoint2D(new Point2D.Double(10, 50));
        polygon.addPoint2D(new Point2D.Double(10, -10));
        System.out.println(polygon);
        System.out.println("");

        Point2D pt = new Point2D.Double(0, 0);
        System.out.println("polygon contains " + pt + ": " + polygon.containsPoint2D(pt));
        pt = new Point2D.Double(100, 100);
        System.out.println("polygon contains " + pt + ": " + polygon.containsPoint2D(pt));
        System.out.println("");

        Rectangle2D boundingBox = new Rectangle2D.Double(0, 0, 100, 100);  // defines lower-left corner and width/height
        System.out.println(boundingBox);
        Polygon2D myClip = (Polygon2D)polygon.clip(boundingBox);
        System.out.println(myClip);

    }

}