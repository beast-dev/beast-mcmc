/*
 * CalendarDateTest.java
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

package test.dr.evoxml;

import dr.evolution.util.Date;
import dr.evolution.util.Units;
import dr.evoxml.DateParser;
import dr.xml.XMLParseException;
import junit.framework.TestCase;

/**
 * Dates given as a value alone: calendar dates, months, years and decimal years.
 */
public class CalendarDateTest extends TestCase {

    public void testCalendarDate() throws Exception {
        // January 1st is the start of the year, in a leap year or not
        assertDate(2020.0, 0.0, "2020-01-01");
        assertDate(2021.0, 0.0, "2021-01-01");
        // the fraction of the year that has passed: 2020 has 366 days, 2021 has 365
        assertDate(2020.0 + 82.0 / 366.0, 0.0, "2020-03-23");
        assertDate(2021.0 + 81.0 / 365.0, 0.0, "2021-03-23");
        assertDate(2020.0 + 365.0 / 366.0, 0.0, "2020-12-31");
        assertDate(2000.0 + 59.0 / 366.0, 0.0, "2000-02-29");
        assertFalse(DateParser.parseCalendarDate("2020-03-23").isUncertaintyCentred());
    }

    public void testMonth() throws Exception {
        // the middle of the month, with the length of the month over 365 as the uncertainty
        assertDate(2021.0 + (0.0 + 31.0) / 2.0 / 365.0, 31.0 / 365.0, "2021-01");
        assertDate(2021.0 + (31.0 + 59.0) / 2.0 / 365.0, 28.0 / 365.0, "2021-02");
        assertDate(2020.0 + (31.0 + 60.0) / 2.0 / 366.0, 29.0 / 365.0, "2020-02");
        assertDate(2021.0 + (334.0 + 365.0) / 2.0 / 365.0, 31.0 / 365.0, "2021-12");
        assertTrue(DateParser.parseCalendarDate("2021-02").isUncertaintyCentred());
    }

    public void testYear() throws Exception {
        assertDate(2020.5, 1.0, "2020");
        assertDate(1918.5, 1.0, "1918");
        assertTrue(DateParser.parseCalendarDate("2020").isUncertaintyCentred());
    }

    public void testDecimalYear() throws Exception {
        assertDate(2020.0, 0.0, "2020.0");
        assertDate(2020.225, 0.0, "2020.225");
        assertDate(2020.0, 0.0, "2020.");
        assertDate(2003.98, 0.0, " 2003.98 ");
    }

    public void testOtherValuesAreLeftAlone() throws Exception {
        // not four digits, so these are read as they always were
        assertNull(DateParser.parseCalendarDate("10.0"));
        assertNull(DateParser.parseCalendarDate("300"));
        assertNull(DateParser.parseCalendarDate("12345"));
        assertNull(DateParser.parseCalendarDate("23/03/2020"));
        assertNull(DateParser.parseCalendarDate("-2020"));
        assertNull(DateParser.parseCalendarDate("2020-3-23"));
    }

    public void testImpossibleDates() {
        assertRejected("2021-02-29");
        assertRejected("2021-02-30");
        assertRejected("2021-13-01");
        assertRejected("2021-00-10");
        assertRejected("2021-13");
        assertRejected("2021-00");
    }

    private static void assertDate(double decimalYear, double uncertainty, String value) throws XMLParseException {
        Date date = DateParser.parseCalendarDate(value);
        assertNotNull(value, date);
        assertEquals(value, decimalYear, date.getTimeValue(), 1E-9);
        assertEquals(value, uncertainty, date.getUncertainty(), 1E-12);
        assertEquals(value, Units.Type.YEARS, date.getUnits());
        assertFalse(value, date.isBackwards());
    }

    private static void assertRejected(String value) {
        try {
            DateParser.parseCalendarDate(value);
            fail(value + " should be rejected");
        } catch (XMLParseException expected) {
        }
    }
}
