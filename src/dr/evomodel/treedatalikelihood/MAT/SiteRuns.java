/*
 * SiteRuns.java
 *
 * Copyright © 2026, the BEAST Development Team.
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
 */

package dr.evomodel.treedatalikelihood.MAT;

import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.Arrays;

/**
 * Sets of patterns held as runs: [start, end) pairs in pattern order, disjoint and not touching.
 */
final class SiteRuns {

    private static final long SIGN_BITS = 0x8080808080808080L;

    private SiteRuns() {
    }

    /**
     * @return the runs of the patterns at which the states are ambiguous (negative)
     */
    static int[] ambiguous(byte[] states, int patternCount) {
        // eight states at a time: a word has an ambiguous (negative) state iff one of its sign bits is set
        final ByteBuffer words = ByteBuffer.wrap(states).order(ByteOrder.nativeOrder());
        final int wordEnd = patternCount - 7;
        int[] runs = new int[8];
        int count = 0;
        int p = 0;
        while (p < patternCount) {
            if (p < wordEnd && (words.getLong(p) & SIGN_BITS) == 0L) {
                p += 8;
                continue;
            }
            if (states[p] < 0) {
                int start = p;
                while (p < patternCount && states[p] < 0) {
                    if (p < wordEnd && (words.getLong(p) & SIGN_BITS) == SIGN_BITS) {
                        p += 8;
                    } else {
                        p++;
                    }
                }
                if (count + 2 > runs.length) {
                    runs = Arrays.copyOf(runs, 2 * runs.length);
                }
                runs[count++] = start;
                runs[count++] = p;
            } else {
                p++;
            }
        }
        return Arrays.copyOf(runs, count);
    }

    /**
     * @return the intersection of two sets of runs
     */
    static int[] intersect(int[] a, int[] b) {
        int[] result = new int[a.length + b.length];
        int count = 0;
        int i = 0;
        int j = 0;
        while (i < a.length && j < b.length) {
            int start = Math.max(a[i], b[j]);
            int end = Math.min(a[i + 1], b[j + 1]);
            if (start < end) {
                result[count++] = start;
                result[count++] = end;
            }
            if (a[i + 1] < b[j + 1]) {
                i += 2;
            } else {
                j += 2;
            }
        }
        return Arrays.copyOf(result, count);
    }

    /**
     * @return the runs of the patterns (of patternCount) not in the given runs
     */
    static int[] complement(int[] runs, int patternCount) {
        int[] result = new int[runs.length + 2];
        int count = 0;
        int position = 0;
        for (int i = 0; i < runs.length; i += 2) {
            if (runs[i] > position) {
                result[count++] = position;
                result[count++] = runs[i];
            }
            position = runs[i + 1];
        }
        if (position < patternCount) {
            result[count++] = position;
            result[count++] = patternCount;
        }
        return Arrays.copyOf(result, count);
    }

    /**
     * @return the patterns in a and not in b
     */
    static int[] subtract(int[] a, int[] b, int patternCount) {
        return intersect(a, complement(b, patternCount));
    }
}
