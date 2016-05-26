/*
 * KMLCoordinates.java
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

import java.util.Arrays;
import java.util.StringTokenizer;

/**
 * @author Marc A. Suchard
 */

public class KMLCoordinates {

    public static final String COORDINATES = "coordinates";
    public static final String FORMAT = "%7.5f";
    public static final String SEPARATOR = ",";
    public static final String NEWLINE = "\n";
    public static final String POINT_SEPARATORS = NEWLINE + " ";

    public KMLCoordinates(double[] x, double[] y) {
        this(x, y, 0.0);
    }

    public KMLCoordinates(double[] x, double[] y, double z) {
        this.x = x;
        this.y = y;

        if (x.length != y.length)
            throw new RuntimeException("Cannot create coordinate system with unbalanced entries");

        this.z = new double[x.length];
        Arrays.fill(this.z, z);
    }

    public KMLCoordinates(double[] x, double[] y, double[] z) {

        this.x = x;
        this.y = y;
        this.z = z;

        if (x.length != y.length && x.length != z.length)
            throw new RuntimeException("Cannot create coordinate system with unbalanced entries");

        length = x.length;
    }

    public void switchXY() {
        double[] tempX = x;
        x = y;
        y = tempX;

    }

    public Element toXML() {
        Element thisElement = new Element(COORDINATES);
        StringBuffer bf = new StringBuffer();
        bf.append(POINT_SEPARATORS);
        for (int i = 0; i < x.length; i++) {
            bf.append(String.format(FORMAT, x[i])).append(SEPARATOR);
            bf.append(String.format(FORMAT, y[i])).append(SEPARATOR);
            bf.append(String.format(FORMAT, z[i])).append(POINT_SEPARATORS);
        }
        thisElement.addContent(bf.toString());
        return thisElement;
    }

    public static XMLObjectParser COORDINATESPARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return COORDINATES;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            String value = (String) xo.getChild(0);

            StringTokenizer st1 = new StringTokenizer(value, POINT_SEPARATORS);
            int count = st1.countTokens();

            double[] x = new double[count];
            double[] y = new double[count];
            double[] z = new double[count];

            for (int i = 0; i < count; i++) {
                String line = st1.nextToken();
                StringTokenizer st2 = new StringTokenizer(line, SEPARATOR);
                if (st2.countTokens() != 3)
                    throw new XMLParseException("All KML coordinates must contain (X,Y,Z) values.  Three dimensions not found in element '" + line + "'");
                x[i] = Double.valueOf(st2.nextToken());
                y[i] = Double.valueOf(st2.nextToken());
                z[i] = Double.valueOf(st2.nextToken());
            }

            return new KMLCoordinates(x, y, z);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element represents a set of (X,Y,Z) coordinates in KML format";
        }

        public Class getReturnType() {
            return KMLCoordinates.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
//            AttributeRule.newDoubleRule(TIME,true),
                new ElementRule(String.class)
        };
    };

    public double[] x;
    public double[] y;
    public double[] z;

    public int length;

}
