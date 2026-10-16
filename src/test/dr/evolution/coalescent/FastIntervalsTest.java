/*
 * FastIntervalsTest.java
 *
 * Copyright © 2002-2026 the BEAST Development Team
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

package test.dr.evolution.coalescent;

import dr.evolution.coalescent.FastIntervals;
import dr.evolution.coalescent.IntervalType;
import junit.framework.TestCase;

import java.util.Arrays;
import java.util.Random;

/**
 * Checks that FastIntervals, which sorts event times incrementally when they are added with node
 * numbers, gives exactly the same intervals as sorting them in full every time.
 */
public class FastIntervalsTest extends TestCase {

    private static final int TIPS = 60;
    private static final int INTERNAL = TIPS - 1;

    private final Random random = new Random(1234);

    public FastIntervalsTest(String name) {
        super(name);
    }

    /**
     * Random updates of a few, some or many times, with store/restore cycles (copyIntervals into a
     * second instance, which is swapped back in on rejection), as TreeIntervals does.
     */
    public void testIncrementalSortingMatchesFullSort() {
        double[] tips = new double[TIPS];
        for (int i = 0; i < TIPS; i++) {
            // several tips share a time, to include ties
            tips[i] = random.nextInt(10) * 0.5;
        }
        double[] internal = new double[INTERNAL];
        for (int i = 0; i < INTERNAL; i++) {
            internal[i] = 5.0 + random.nextDouble() * 20.0;
        }

        FastIntervals current = new FastIntervals(TIPS, INTERNAL);
        FastIntervals stored = new FastIntervals(TIPS, INTERNAL);
        addEvents(current, tips, internal, true);
        assertSameAsFullSort(current, tips, internal);

        for (int iteration = 0; iteration < 5000; iteration++) {
            // store
            stored.copyIntervals(current);
            double[] storedTips = tips.clone();
            double[] storedInternal = internal.clone();

            // propose: change none, one, a few or many times (occasionally tip times, or a tie)
            int changes = random.nextInt(10) == 0 ? random.nextInt(INTERNAL) : random.nextInt(4);
            for (int c = 0; c < changes; c++) {
                int i = random.nextInt(INTERNAL);
                internal[i] = random.nextInt(5) == 0 ? internal[random.nextInt(INTERNAL)] : 5.0 + random.nextDouble() * 20.0;
            }
            if (random.nextInt(20) == 0) {
                tips[random.nextInt(TIPS)] = random.nextInt(10) * 0.5;
            }
            addEvents(current, tips, internal, true);
            assertSameAsFullSort(current, tips, internal);

            if (random.nextBoolean()) {
                // reject: restore by swapping, then the next evaluation uses the restored instance
                FastIntervals tmp = current;
                current = stored;
                stored = tmp;
                tips = storedTips;
                internal = storedInternal;
                addEvents(current, tips, internal, true);
                assertSameAsFullSort(current, tips, internal);
            }
        }
    }

    /**
     * Events added without node numbers (in any order) are fully sorted.
     */
    public void testEventsWithoutNodeNumbers() {
        double[] tips = new double[TIPS];
        double[] internal = new double[INTERNAL];
        FastIntervals intervals = new FastIntervals(TIPS, INTERNAL);
        for (int iteration = 0; iteration < 200; iteration++) {
            for (int i = 0; i < TIPS; i++) {
                tips[i] = random.nextInt(10) * 0.5;
            }
            for (int i = 0; i < INTERNAL; i++) {
                internal[i] = 5.0 + random.nextDouble() * 20.0;
            }
            // alternate between the two ways of adding events
            addEvents(intervals, tips, internal, iteration % 3 != 0);
            assertSameAsFullSort(intervals, tips, internal);
        }
    }

    private static void addEvents(FastIntervals intervals, double[] tips, double[] internal, boolean withNodeNumbers) {
        intervals.resetEvents();
        for (int i = 0; i < tips.length; i++) {
            if (withNodeNumbers) {
                intervals.addSampleEvent(tips[i], i);
            } else {
                intervals.addSampleEvent(tips[i]);
            }
        }
        for (int i = 0; i < internal.length; i++) {
            if (withNodeNumbers) {
                intervals.addCoalescentEvent(internal[i], tips.length + i);
            } else {
                intervals.addCoalescentEvent(internal[i]);
            }
        }
    }

    /**
     * Compares with the intervals given by fully sorting the times (the original algorithm).
     */
    private static void assertSameAsFullSort(FastIntervals intervals, double[] tips, double[] internal) {
        double[] sampleTimes = tips.clone();
        double[] coalescentTimes = internal.clone();
        Arrays.sort(sampleTimes);
        Arrays.sort(coalescentTimes);

        int count = tips.length + internal.length - 1;
        double[] expectedIntervals = new double[count];
        double[] expectedTimes = new double[count];
        IntervalType[] expectedTypes = new IntervalType[count];
        int[] expectedLineages = new int[count];
        double lastTime = sampleTimes[0];
        int s = 1, c = 0, i = 0, lineages = 1;
        while (s < sampleTimes.length || c < coalescentTimes.length) {
            if (s < sampleTimes.length && sampleTimes[s] <= coalescentTimes[c]) {
                expectedIntervals[i] = sampleTimes[s] - lastTime;
                expectedTimes[i] = sampleTimes[s];
                expectedTypes[i] = IntervalType.SAMPLE;
                expectedLineages[i] = lineages;
                lastTime = sampleTimes[s];
                lineages++;
                s++;
            } else {
                expectedIntervals[i] = coalescentTimes[c] - lastTime;
                expectedTimes[i] = coalescentTimes[c];
                expectedTypes[i] = IntervalType.COALESCENT;
                expectedLineages[i] = lineages;
                lastTime = coalescentTimes[c];
                lineages--;
                c++;
            }
            i++;
        }

        assertEquals(count, intervals.getIntervalCount());
        for (int k = 0; k < count; k++) {
            // exact equality: the incremental sort must give the same values in the same order
            assertEquals("interval " + k, expectedIntervals[k], intervals.getInterval(k), 0.0);
            assertEquals("time " + k, expectedTimes[k], intervals.getIntervalTime(k), 0.0);
            assertEquals("type " + k, expectedTypes[k], intervals.getIntervalType(k));
            assertEquals("lineages " + k, expectedLineages[k], intervals.getLineageCount(k));
        }
        assertEquals(sampleTimes[0], intervals.getStartTime(), 0.0);
    }
}
