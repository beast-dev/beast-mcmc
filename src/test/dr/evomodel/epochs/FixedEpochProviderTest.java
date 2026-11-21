/*
 * FixedEpochProviderTest.java
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

package test.dr.evomodel.epochs;

import dr.evolution.util.Units;
import dr.evomodel.epochs.EpochProvider;
import dr.evomodel.epochs.FixedEpochProvider;
import dr.evomodelxml.epochs.FixedEpochProviderParser;
import junit.framework.TestCase;

public class FixedEpochProviderTest extends TestCase {

    public void testEqualDurations() {
        EpochProvider epochs = FixedEpochProvider.equalDurations(5, 0.5);
        assertEquals(5, epochs.getEpochCount());
        assertEquals(4, epochs.getTransitionTimes().length);
        for (int i = 0; i < 4; i++) {
            assertEquals((i + 1) * 0.5, epochs.getTransitionTimes()[i], 0.0);
        }
    }

    public void testOneEpoch() {
        EpochProvider epochs = FixedEpochProvider.equalDurations(1, 2.0);
        assertEquals(1, epochs.getEpochCount());
        assertEquals(0, epochs.getTransitionTimes().length);
        assertEquals(0, epochs.getTransitionTimeParameter().getDimension());
        assertEquals(0, epochs.getEpochIndex(0.0));
        assertEquals(0, epochs.getEpochIndex(1000.0));
    }

    public void testGivenTransitionTimes() {
        double[] times = {0.1, 0.25, 1.0, 7.5};
        EpochProvider epochs = new FixedEpochProvider(times);
        assertEquals(5, epochs.getEpochCount());
        for (int i = 0; i < times.length; i++) {
            assertEquals(times[i], epochs.getTransitionTimes()[i], 0.0);
        }
        // the provider keeps its own copy
        times[0] = 99.0;
        assertEquals(0.1, epochs.getTransitionTimes()[0], 0.0);
    }

    public void testParameterMatchesTimes() {
        EpochProvider epochs = new FixedEpochProvider(new double[]{0.1, 0.25, 1.0});
        assertEquals(3, epochs.getTransitionTimeParameter().getDimension());
        for (int i = 0; i < 3; i++) {
            assertEquals(epochs.getTransitionTimes()[i], epochs.getTransitionTimeParameter().getParameterValue(i), 0.0);
        }
    }

    public void testEpochIndex() {
        EpochProvider epochs = new FixedEpochProvider(new double[]{1.0, 2.0, 4.0});
        assertEquals(0, epochs.getEpochIndex(0.0));
        assertEquals(0, epochs.getEpochIndex(0.999));
        // a transition time is in the older epoch
        assertEquals(1, epochs.getEpochIndex(1.0));
        assertEquals(1, epochs.getEpochIndex(1.5));
        assertEquals(2, epochs.getEpochIndex(2.0));
        assertEquals(2, epochs.getEpochIndex(3.999));
        assertEquals(3, epochs.getEpochIndex(4.0));
        assertEquals(3, epochs.getEpochIndex(1.0E6));
    }

    public void testBadTransitionTimes() {
        assertRejected(new double[]{1.0, 0.5});
        assertRejected(new double[]{1.0, 1.0});
        assertRejected(new double[]{0.0, 1.0});
        assertRejected(new double[]{-1.0, 1.0});
        assertRejected(new double[]{1.0, Double.NaN});
        assertRejected(new double[]{1.0, Double.POSITIVE_INFINITY});
    }

    public void testBadEqualDurations() {
        try {
            FixedEpochProvider.equalDurations(0, 1.0);
            fail("no epochs should be rejected");
        } catch (IllegalArgumentException expected) {
        }
        try {
            FixedEpochProvider.equalDurations(3, 0.0);
            fail("a duration of zero should be rejected");
        } catch (IllegalArgumentException expected) {
        }
    }

    public void testUnitConversion() {
        double[] one = {1.0};
        assertEquals(1.0 / 12.0, convert(one, Units.Type.MONTHS, Units.Type.YEARS), 1E-15);
        assertEquals(7.0 / 365.25, convert(one, Units.Type.WEEKS, Units.Type.YEARS), 1E-15);
        assertEquals(1.0 / 365.25, convert(one, Units.Type.DAYS, Units.Type.YEARS), 1E-15);
        assertEquals(1.0, convert(one, Units.Type.YEARS, Units.Type.YEARS), 0.0);
        assertEquals(365.25, convert(one, Units.Type.YEARS, Units.Type.DAYS), 1E-12);
        assertEquals(7.0, convert(one, Units.Type.WEEKS, Units.Type.DAYS), 0.0);
        assertEquals(365.25 / 7.0, convert(one, Units.Type.YEARS, Units.Type.WEEKS), 1E-12);
        assertEquals(365.25 / 12.0 / 7.0, convert(one, Units.Type.MONTHS, Units.Type.WEEKS), 1E-12);
        assertEquals(1.0 / 7.0, convert(one, Units.Type.DAYS, Units.Type.WEEKS), 1E-15);

        // 12 epochs of a month in a tree in years
        double[] months = FixedEpochProvider.equalDurationTransitionTimes(12, 1.0);
        double[] years = FixedEpochProviderParser.convert(months, Units.Type.MONTHS, Units.Type.YEARS);
        assertEquals(11, years.length);
        assertEquals(0.5, years[5], 1E-15);
        assertEquals(11.0 / 12.0, years[10], 1E-15);
    }

    private static double convert(double[] times, Units.Type from, Units.Type to) {
        return FixedEpochProviderParser.convert(times, from, to)[0];
    }

    private static void assertRejected(double[] transitionTimes) {
        try {
            new FixedEpochProvider(transitionTimes);
            fail("should be rejected: " + java.util.Arrays.toString(transitionTimes));
        } catch (IllegalArgumentException expected) {
        }
    }
}
