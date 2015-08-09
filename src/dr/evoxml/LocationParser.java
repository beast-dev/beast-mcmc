/*
 * LocationParser.java
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

package dr.evoxml;

import dr.xml.*;
import dr.evolution.util.Location;

/**
 * @author Andrew Rambaut
 *
 * @version $Id: DateParser.java,v 1.2 2005/05/24 20:25:59 rambaut Exp $
 */
public class LocationParser extends AbstractXMLObjectParser {

    public static final String DESCRIPTION = "description";
    public static final String LONGITUDE = "longitude";
    public static final String LATITUDE = "latitude";

    public String getParserName() { return Location.LOCATION; }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        if (xo.getChildCount() > 0) {
            throw new XMLParseException("No child elements allowed in location element.");
        }

        String description = xo.getAttribute(DESCRIPTION, "");

        double longitude = parseLongLat(xo.getAttribute(LONGITUDE, ""));
        double latitude = parseLongLat(xo.getAttribute(LATITUDE, ""));

        return Location.newLocation(xo.getId(), description, longitude, latitude);
    }

    private double parseLongLat(final String value) throws XMLParseException {
        double d = 0.0;

        if (value != null && value.length() > 0) {
            try {
                d = Double.parseDouble(value);
            } catch (NumberFormatException nfe) {
                // @todo - parse degrees minutes and seconds
            }
        }

        return d;
    }

    public String getParserDescription() {
        return "Specifies a location with an optional longitude and latitude";
    }

    public XMLSyntaxRule[] getSyntaxRules() { return rules; }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
            AttributeRule.newStringRule(DESCRIPTION, true,
                    "A description of this location"),
            AttributeRule.newStringRule(LONGITUDE, true,
                    "The longitude in degrees, minutes, seconds or decimal degrees"),
            AttributeRule.newStringRule(LATITUDE, true,
                    "The latitude in degrees, minutes, seconds or decimal degrees"),
    };

    public Class getReturnType() { return Location.class; }
}