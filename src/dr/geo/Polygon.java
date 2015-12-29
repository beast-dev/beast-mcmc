/*
 * Polygon.java
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

import org.jdom.Element;

import java.util.StringTokenizer;
import java.util.List;

import dr.xml.XMLParseException;

/**
 * Created by IntelliJ IDEA.
 * User: phil
 * Date: Feb 19, 2009
 * Time: 5:15:25 PM
 * To change this template use File | Settings | File Templates.
 */
public class Polygon {

    public static final String POLYGON = "polygon";
    public static final String CLOSED = "closed";
    public static final String COORDINATES = "coordinates";
    // includes a density value
    //public static final String DENSITY = "density";

    public Polygon(Element e) {

        boolean closed = Boolean.valueOf(e.getAttributeValue(CLOSED));

        List children = e.getChildren();
        for (int a = 0; a < children.size(); a++) {
            Element childElement = (Element)children.get(a);
            if (childElement.getName().equals(COORDINATES)) {

            String value = childElement.getTextTrim();
            StringTokenizer st1 = new StringTokenizer(value, "\n");
            int count = st1.countTokens();  //System.out.println(count);

            x = new double[count];
            y = new double[count];
            z = new double[count];
            // assumes closed!
            length = x.length - 1;

            for (int i = 0; i < count; i++) {
                String line = st1.nextToken();
                StringTokenizer st2 = new StringTokenizer(line, ",");
                if (st2.countTokens() != 3)
                    throw new IllegalArgumentException("All KML coordinates must contain (X,Y,Z) values.  Three dimensions not found in element '" + line + "'");
                x[i] = Double.valueOf(st2.nextToken());
                y[i] = Double.valueOf(st2.nextToken());
                z[i] = Double.valueOf(st2.nextToken());
                }
            //gets the density value
            //String densityValue = e.getChild(DENSITY).getTextTrim();
            //d = Double.valueOf(densityValue);
            }
        }
        planarInZ = isPlanarInZ();
    }

    protected double[] x;
    protected double[] y;
    protected double[] z;
    protected double[] t;
    // declares dernsity double
    //protected double d;

    protected int length;
    private boolean planarInZ;

    //returns density value
    //public double getDensity(){
        //return d;
    //}

    private boolean isPlanarInZ() {
        double z0 = z[0];
        boolean planar = true;
        for (int i = 0; i < length; i++) {
            if (z0 != z[i]) {
                planar = false;
                break;
            }
        }
        return planar;
    }

    public boolean contains2DPoint(double inX, double inY) {

//        System.err.print("contains (" + inX + "," + inY + ") = ");

        if (!planarInZ)
            throw new RuntimeException("Only 2D polygons are currently implemented");
        boolean contains = false;

        // Take a horizontal ray from (inX,inY) to the right.
        // If ray across the polygon edges an odd # of times, the point is inside.
        for (int i = 0, j = length - 1; i < length; j = i++) {

             if ((((y[i] <= inY) && (inY < y[j])) ||
                    ((y[j] <= inY) && (inY < y[i]))) &&
                    (inX < (x[j] - x[i]) * (inY - y[i]) / (y[j] - y[i]) + x[i]))
                contains = !contains;
        }

//        System.err.println(contains);
        return contains;
    }

}
