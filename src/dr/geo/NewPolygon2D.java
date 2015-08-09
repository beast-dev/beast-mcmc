/*
 * NewPolygon2D.java
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

import dr.xml.*;
import org.jdom.Element;

import java.awt.geom.*;
import java.util.List;
import java.util.StringTokenizer;

/**
 * @author Marc A. Suchard
 * @author Philippe Lemey
 */

public class NewPolygon2D {

    public static final String POLYGON = "polygon";
    public static final String CLOSED = "closed";
    public static final String FILL_VALUE = "fillValue";

    public NewPolygon2D(GeneralPath path) {
        this.path = path;
    }

    public NewPolygon2D(Element e) {

        List<Element> children = e.getChildren();
//        for (int a = 0; a < children.size(); a++) {
        for (Element childElement : children) {
//            Element childElement = (Element) children.get(a);
            if (childElement.getName().equals(KMLCoordinates.COORDINATES)) {

                String value = childElement.getTextTrim();
                StringTokenizer st1 = new StringTokenizer(value, "\n");
                int count = st1.countTokens();  //System.out.println(count);

//                point2Ds = new LinkedList<Point2D>();
                this.path = new GeneralPath();

                for (int i = 0; i < count; i++) {
                    String line = st1.nextToken();
                    StringTokenizer st2 = new StringTokenizer(line, ",");
                    if (st2.countTokens() != 3)
                        throw new IllegalArgumentException("All KML coordinates must contain (X,Y,Z) values.  Three dimensions not found in element '" + line + "'");
                    final float x = Float.valueOf(st2.nextToken());
                    final float y = Float.valueOf(st2.nextToken());
                    if (i == 0) {
                        path.moveTo(x, y);
                    } else {
                        path.lineTo(x, y);
                    }
//                    point2Ds.add(new Point2D.Double(x, y));
                }
//                length = point2Ds.size() - 1;
                break;

            }
        }
    }

    public NewPolygon2D() {
        path = new GeneralPath();
    }

    public void moveTo(Point2D pt) {
        path.moveTo((float) pt.getX(), (float) pt.getY());
    }

    public void lineTo(Point2D pt) {
        path.lineTo((float) pt.getX(), (float) pt.getY());
    }

    public boolean contains(Point2D pt) {
        return path.contains(pt);
    }

    public void closePath() {
        path.closePath();
    }

    public void setFillValue(double value) {
        fillValue = value;
    }

    public double getFillValue() {
        return fillValue;
    }

    public NewPolygon2D clip(Rectangle2D boundingBox) {
        Area thisArea = new Area(path);
        thisArea.intersect(new Area(boundingBox));
        PathIterator iterator = thisArea.getPathIterator(null);
        double[] v = new double[2];
        while (!iterator.isDone()) {
            int type = iterator.currentSegment(v);
            System.err.println(":" + v[0] + v[1] + "\n");
            iterator.next();
        }
        System.exit(-1);

        GeneralPath path = new GeneralPath(thisArea);
        path.closePath();
        NewPolygon2D newPolygon = new NewPolygon2D(path);
        newPolygon.setFillValue(this.getFillValue());
        return newPolygon;
    }

    public String toString() {
        StringBuffer sb = new StringBuffer();
        sb.append("polygon(");
        sb.append(fillValue);
        sb.append(")[\n");
        PathIterator iterator = path.getPathIterator(null);
        float[] values = new float[2];
        while (!iterator.isDone()) {
            int type = iterator.currentSegment(values);
            Point2D pt = new Point2D.Double(values[0], values[1]);
            sb.append("\t");
            sb.append(pt);
            sb.append("\n");
            iterator.next();
        }
//        sb.append(iterator);
//        for(Point2D pt : point2Ds) {
//            sb.append("\t");
//            sb.append(pt);
//            sb.append("\n");
//        }
//        sb.append(path.toString());
        sb.append("]");
        return sb.toString();
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return POLYGON;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            KMLCoordinates coordinates = (KMLCoordinates) xo.getChild(KMLCoordinates.class);
            boolean closed = xo.getAttribute(CLOSED, false);

            if ((!closed && coordinates.length < 3) ||
                    (closed && coordinates.length < 4))
                throw new XMLParseException("Insufficient points in polygon '" + xo.getId() + "' to define a polygon in 2D");

            NewPolygon2D polygon = new NewPolygon2D();
            polygon.moveTo(new Point2D.Double(coordinates.x[0], coordinates.y[0]));
            int length = coordinates.length;
            if (closed)
                length--;
            for (int i = 1; i < length; i++)
                polygon.lineTo(new Point2D.Double(coordinates.x[i], coordinates.y[i]));
            polygon.lineTo(new Point2D.Double(coordinates.x[0], coordinates.y[0]));
//            polygon.closePath();

            polygon.setFillValue(xo.getAttribute(FILL_VALUE, 0.0));

            return polygon;
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element represents a polygon.";
        }

        public Class getReturnType() {
            return Polygon2D.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(KMLCoordinates.class),
                AttributeRule.newBooleanRule(CLOSED, true),
                AttributeRule.newDoubleRule(FILL_VALUE, true),
        };
    };

    public static void main(String[] args) {
        NewPolygon2D polygon = new NewPolygon2D();
        polygon.moveTo(new Point2D.Double(-10, -10));
        polygon.lineTo(new Point2D.Double(-10, 50));
        polygon.lineTo(new Point2D.Double(10, 50));
        polygon.lineTo(new Point2D.Double(10, -10));
        polygon.lineTo(new Point2D.Double(-10, -10));
//        polygon.closePath();
        System.out.println(polygon);
        System.out.println("");
//        System.exit(-1);

        Point2D pt = new Point2D.Double(0, 0);
        System.out.println("polygon contains " + pt + ": " + polygon.contains(pt));
        pt = new Point2D.Double(100, 100);
        System.out.println("polygon contains " + pt + ": " + polygon.contains(pt));
        System.out.println("");

        Rectangle2D boundingBox = new Rectangle2D.Double(0, 0, 100, 100);  // defines lower-left corner and width/height
        System.out.println(boundingBox);
        NewPolygon2D myClip = polygon.clip(boundingBox);
        System.out.println(myClip);

    }

    private double fillValue;
    private GeneralPath path;
}