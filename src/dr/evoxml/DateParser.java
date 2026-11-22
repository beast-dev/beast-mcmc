/*
 * DateParser.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.evoxml;

import dr.evolution.util.Date;
import dr.evolution.util.Units;
import dr.evoxml.util.DateUnitsType;
import dr.evoxml.util.XMLUnits;
import dr.xml.*;

import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.TimeZone;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
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

        if (xo.getChildCount() > 0) {
            throw new XMLParseException("No child elements allowed in date element.");
        }

        if (xo.hasAttribute(VALUE) && !xo.hasAttribute(DIRECTION) && !xo.hasAttribute(UNITS) &&
                !xo.hasAttribute(UNCERTAINTY) && !xo.hasAttribute(PRECISION) && !xo.hasAttribute(ORIGIN)) {
            // a value with nothing else can be a calendar date or a decimal year
            Date calendarDate = parseCalendarDate(xo.getStringAttribute(VALUE));
            if (calendarDate != null) {
                return calendarDate;
            }
        }

        double value = 0.0;
        java.util.Date dateValue = null;

        if (xo.hasAttribute(VALUE)) {
            try {
                value = xo.getDoubleAttribute(VALUE);
            } catch (XMLParseException e) {
                // not a number, so it has to be a calendar date
                dateValue = parseCalendarDay(xo.getStringAttribute(VALUE), VALUE);
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

            java.util.Date origin = parseCalendarDay(xo.getStringAttribute(ORIGIN), ORIGIN);

            if (dateValue != null) {
                date = new Date(dateValue, units, origin);
            } else {
                date = new Date(value, units, backwards, origin);
            }

        } else {

            if (dateValue != null && units == Units.Type.YEARS) {
                // a calendar date in years is its decimal year, as it is when given with no other attributes
                date = new Date(new Date(dateValue).getTimeValue(), units, false);
            } else if (dateValue != null) {
                // in other units, the time since the default origin (1st Jan 1970)
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

    /**
     * Reads a date given as a value alone, when it is in one of these forms:
     * <ul>
     * <li>yyyy-mm-dd, a calendar date: the start of that day;</li>
     * <li>yyyy-mm, a month: the middle of the month, with an uncertainty of the length of the month;</li>
     * <li>yyyy, a year: the middle of the year, with an uncertainty of a year;</li>
     * <li>yyyy.x, a decimal year (yyyy.0 is the start of January 1st).</li>
     * </ul>
     * The date is in years, running forwards. Where there is an uncertainty the date is the middle of it.
     *
     * @param value the value of the date
     * @return the date, or null if the value is in none of these forms
     * @throws XMLParseException if the value is in one of these forms but is not a date (2001-02-30)
     */
    public static Date parseCalendarDate(String value) throws XMLParseException {
        value = value.trim();

        if (value.matches("\\d\\d\\d\\d-\\d\\d-\\d\\d")) {
            return createDate(calendarDateToDecimalYear(value), 0.0);
        }

        if (value.matches("\\d\\d\\d\\d-\\d\\d")) {
            final int year = Integer.parseInt(value.substring(0, 4));
            final int month = Integer.parseInt(value.substring(5, 7));
            if (month < 1 || month > 12) {
                throw new XMLParseException(VALUE + "=\"" + value + "\" is not a month: use YYYY-MM");
            }
            final double start = calendarDateToDecimalYear(value + "-01");
            final double end = month == 12 ? year + 1.0 :
                    calendarDateToDecimalYear(String.format("%04d-%02d-01", year, month + 1));
            final Calendar calendar = Calendar.getInstance(TimeZone.getTimeZone("GMT"));
            calendar.clear();
            calendar.set(year, month - 1, 1);
            final int daysInMonth = calendar.getActualMaximum(Calendar.DAY_OF_MONTH);
            return createDate((start + end) / 2.0, daysInMonth / 365.0);
        }

        if (value.matches("\\d\\d\\d\\d")) {
            return createDate(Integer.parseInt(value) + 0.5, 1.0);
        }

        if (value.matches("\\d\\d\\d\\d\\.\\d*")) {
            return createDate(Double.parseDouble(value), 0.0);
        }

        return null;
    }

    /**
     * @param value a calendar date, yyyy-mm-dd
     * @return the start of that day as a decimal year: the year plus the fraction of that year that has
     * passed, so January 1st is yyyy.0 whether or not it is a leap year. This is the conversion BEAUti
     * uses for the dates of tips.
     */
    public static double calendarDateToDecimalYear(String value) throws XMLParseException {
        return new Date(parseCalendarDay(value, VALUE)).getTimeValue();
    }

    /**
     * Reads a calendar date. The only form read is YYYY-MM-DD.
     *
     * @param value     the calendar date
     * @param attribute the name of the attribute it is the value of, for the message if it is not a date
     * @return the start of that day (in GMT, so that the date does not depend on where it is read)
     */
    private static java.util.Date parseCalendarDay(String value, String attribute) throws XMLParseException {
        value = value.trim();
        if (value.matches("\\d\\d\\d\\d-\\d\\d-\\d\\d")) {
            final SimpleDateFormat format = new SimpleDateFormat("yyyy-MM-dd");
            format.setLenient(false);
            format.setTimeZone(TimeZone.getTimeZone("GMT"));
            try {
                return format.parse(value);
            } catch (ParseException pe) {
                // the right form but not a day there has been (2021-02-30)
            }
        }
        throw new XMLParseException(attribute + "=\"" + value + "\" is not a calendar date: use YYYY-MM-DD" +
                (value.matches("\\d+/\\d+/\\d+") ? " (dates as DD/MM/YYYY are no longer read)" : ""));
    }

    private static Date createDate(double decimalYear, double uncertainty) {
        Date date = new Date(decimalYear, Units.Type.YEARS, false);
        if (uncertainty > 0.0) {
            date.setUncertainty(uncertainty, true);
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
                        "<date value=\"300.0\" origin=\"1989-01-01\" units=\"days\" direction=\"forwards\"/>\n" +
                        "\n" +
                        "<!-- calendar dates and decimal years, given as a value alone                 -->\n" +
                        "<date value=\"2020-03-23\"/>\n" +
                        "<date value=\"2020-03\"/>   <!-- some time in March 2020 -->\n" +
                        "<date value=\"2020\"/>      <!-- some time in 2020 -->\n" +
                        "<date value=\"2020.225\"/>\n";
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new StringAttributeRule(VALUE,
                    "The value of this date. With no other attributes it can be a calendar date (YYYY-MM-DD), " +
                            "a month (YYYY-MM), a year (YYYY) or a decimal year (YYYY.x)"),
            new StringAttributeRule(ORIGIN,
                    "The origin of this time scale, which must be a calendar date as YYYY-MM-DD", "2001-01-01", true),
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
