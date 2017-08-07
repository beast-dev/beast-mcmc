/*
 * DateParser.java
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

import dr.evolution.util.Date;
import dr.evolution.util.Units;
import dr.evoxml.util.DateUnitsType;
import dr.evoxml.util.XMLUnits;
import dr.xml.*;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: DateParser.java,v 1.2 2005/05/24 20:25:59 rambaut Exp $
 */
public class DateParser extends AbstractXMLObjectParser {

    public static final String VALUE = "value";
    public static final String UNITS = "units";
    public static final String ORIGIN = "origin";
    public static final String DIRECTION = "direction";

    public static final String FORWARDS = DateUnitsType.FORWARDS.getAttribute(); //"forwards";
    public static final String BACKWARDS = DateUnitsType.BACKWARDS.getAttribute(); //"backwards";

    public static final String YEARS = DateUnitsType.YEARS.getAttribute(); //"units";
    public static final String MONTHS = DateUnitsType.MONTHS.getAttribute(); //"units";
    public static final String DAYS = DateUnitsType.DAYS.getAttribute(); //"days";

    public static final String UNCERTAINTY = "uncertainty";
    public static final String PRECISION = "precision"; // inappropriate older attribute

    public String getParserName() {
        return Date.DATE;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        java.text.DateFormat dateFormat = java.text.DateFormat.getDateInstance(java.text.DateFormat.SHORT, java.util.Locale.UK);
        dateFormat.setLenient(true);

        if (xo.getChildCount() > 0) {
            throw new XMLParseException("No child elements allowed in date element.");
        }

        double value = 0.0;
        java.util.Date dateValue = null;

        if (xo.hasAttribute(VALUE)) {
            try {
                value = xo.getDoubleAttribute(VALUE);
            } catch (XMLParseException e) {
                String dateString = xo.getStringAttribute(VALUE);

                try {

                    dateValue = dateFormat.parse(dateString);

                } catch (Exception ex) {
                    throw new XMLParseException("value=" + dateString + " not recognised as a date, use DD/MM/YYYY");
                }
            }
        } else {
            throw new XMLParseException("Value attribute missing from date element.");
        }

        boolean backwards = false;

        if (xo.hasAttribute(DIRECTION)) {
            String direction = (String) xo.getAttribute(DIRECTION);
            if (direction.equals(BACKWARDS)) {
                backwards = true;
            }
        }

        Units.Type units = XMLUnits.Utils.getUnitsAttr(xo);

        Date date;

        if (xo.hasAttribute(ORIGIN)) {

            String originString = (String) xo.getAttribute(ORIGIN);
            java.util.Date origin;

            try {
                origin = dateFormat.parse(originString);
            } catch (Exception e) {
                throw new XMLParseException("origin=" + originString + " not recognised as a date, use DD/MM/YYYY");
            }

            if (dateValue != null) {
                date = new Date(dateValue, units, origin);
            } else {
                date = new Date(value, units, backwards, origin);
            }

        } else {

            // No origin specified so use default (1st Jan 1970)
            if (dateValue != null) {
                date = new Date(dateValue, units);
            } else {
                date = new Date(value, units, backwards);
            }
        }

        if (xo.hasAttribute(UNCERTAINTY)) {
            double uncertainty = xo.getDoubleAttribute(UNCERTAINTY);
            date.setUncertainty(uncertainty);
        } else if (xo.hasAttribute(PRECISION)) {
            // Precision was the old term for uncertainty which was confusing
            // because it was the reciprocal. Parse this term synonymously for
            // backwards compatibility.
            double uncertainty = xo.getDoubleAttribute(PRECISION);
            date.setUncertainty(uncertainty);
        }


        return date;
    }

    public String getParserDescription() {
        return "Specifies a date on a given timescale";
    }

    public String getExample() {
        return
                "<!-- a date representing 10 years in the past                                 -->\n" +
                        "<date value=\"10.0\" units=\"years\" direction=\"backwards\"/>\n" +
                        "\n" +
                        "<!-- a date representing 300 days after Jan 1st 1989                          -->\n" +
                        "<date value=\"300.0\" origin=\"01/01/89\" units=\"days\" direction=\"forwards\"/>\n";
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new StringAttributeRule(VALUE,
                    "The value of this date"),
            new StringAttributeRule(ORIGIN,
                    "The origin of this time scale, which must be a valid calendar date", "01/01/01", true),
            new StringAttributeRule(UNITS, "The units of the timescale", new String[]{YEARS, MONTHS, DAYS}, true),
            new StringAttributeRule(DIRECTION, "The direction of the timescale", new String[]{FORWARDS, BACKWARDS}, true),
            new XORRule(
                    AttributeRule.newDoubleRule(UNCERTAINTY, false, "The uncertainty to which the date is specified"),
                    AttributeRule.newDoubleRule(PRECISION, false, "The uncertainty to which the date is specified"), true
            )
    };

    public Class getReturnType() {
        return dr.evolution.util.Date.class;
    }
}
