/*
 * Polygon2DFill.java
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

import org.jdom.Element;
import java.awt.geom.Point2D;

/**
 * 2D Polygon that uses its fillValue attribute as inclusion probability.
 * These polygons are typically non-overlapping and a list of such polygons needs its fillValues to sum to one.
 * @author Guy Baele
 */
public class Polygon2DFill extends Polygon2D {

    public Polygon2DFill(double[] x, double[] y, double fillValue) {
        super(x,y);
        this.fillValue = fillValue;
    }

    public Polygon2DFill(Element e, double fillValue) {
        super(e);

        this.fillValue = fillValue;
        //no need to recompute this all the time
        this.logFillValue = Math.log(fillValue);
        //this.hasFillValue = true;
    }

    public double getProbability(Point2D Point2D, boolean outside) {
        boolean contains = containsPoint2D(Point2D);
        if (contains) {
            return this.fillValue;
        } else {
            return 0.0;
        }
    }

    public double getLogProbability(Point2D Point2D, boolean outside) {
        boolean contains = containsPoint2D(Point2D);
        if (contains) {
            return this.logFillValue;
        } else {
            return Double.NEGATIVE_INFINITY;
        }
    }

    @Override
    public void setFillValue(double value) {
        fillValue = value;
        logFillValue = Math.log(value);
        //hasFillValue = true;
    }

    @Override
    public double getFillValue() {
        return fillValue;
    }

    @Override
    public double getLogFillValue() {
        return logFillValue;
    }

    @Override
    public boolean hasFillValue() {
        return true;
    }

    private double fillValue;
    private double logFillValue;
    //private boolean hasFillValue;

    public static void main(String[] args) {

    }

}
