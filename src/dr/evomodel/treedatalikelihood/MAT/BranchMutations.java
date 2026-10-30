/*
 * BranchMutations.java
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

import java.util.Arrays;

/**
 * The events on one branch, held in primitive arrays: the substitutions first, then the
 * ambiguities (sites becoming ambiguous, with a to-state of -1). Instances are immutable, so a
 * branch's events can be cached against the identity of its BranchMutations object.
 *
 * With a single-hit branch history model there is at most one substitution per site. Models with
 * explicit multiple hits store a site's hits contiguously, in order along the branch (from the parent
 * end), so that the hits of a site form a path from the parent's state to the child's. Models with
 * timed hits also store each hit's position as a fraction of the branch length from the parent end.
 */
public final class BranchMutations {

    public static final BranchMutations EMPTY = new BranchMutations(new int[0], new byte[0], new byte[0], 0, null);

    private final int[] patterns;
    private final byte[] fromStates;
    private final byte[] toStates;
    private final int substitutionCount;
    private final double[] times;

    BranchMutations(int[] patterns, byte[] fromStates, byte[] toStates, int substitutionCount) {
        this(patterns, fromStates, toStates, substitutionCount, null);
    }

    BranchMutations(int[] patterns, byte[] fromStates, byte[] toStates, int substitutionCount, double[] times) {
        this.patterns = patterns;
        this.fromStates = fromStates;
        this.toStates = toStates;
        this.substitutionCount = substitutionCount;
        this.times = times;
    }

    /**
     * @return the number of events (substitutions and ambiguities)
     */
    public int size() {
        return patterns.length;
    }

    /**
     * @return the number of substitutions, which are events 0 to getSubstitutionCount() - 1
     */
    public int getSubstitutionCount() {
        return substitutionCount;
    }

    public int getAmbiguityCount() {
        return patterns.length - substitutionCount;
    }

    public int getPattern(int i) {
        return patterns[i];
    }

    public int getFromState(int i) {
        return fromStates[i];
    }

    /**
     * @return the state below the event, or -1 for an ambiguity
     */
    public int getToState(int i) {
        return toStates[i];
    }

    public boolean isAmbiguity(int i) {
        return i >= substitutionCount;
    }

    /**
     * @return true if the substitutions carry positions along the branch
     */
    public boolean hasTimes() {
        return times != null;
    }

    /**
     * @return the position of substitution i as a fraction of the branch length from the parent end
     */
    public double getTime(int i) {
        return times[i];
    }

    /**
     * @return the index after the last substitution at the same site as substitution i (the hits of a
     * site are contiguous)
     */
    public int getSiteEnd(int i) {
        int end = i + 1;
        while (end < substitutionCount && patterns[end] == patterns[i]) {
            end++;
        }
        return end;
    }

    /**
     * @return the events as Mutation objects (a convenience; not for use in tight loops)
     */
    public Mutation[] toMutations(int nodeNumber) {
        Mutation[] mutations = new Mutation[patterns.length];
        for (int i = 0; i < patterns.length; i++) {
            mutations[i] = isAmbiguity(i) ?
                    Mutation.ambiguity(nodeNumber, patterns[i], fromStates[i]) :
                    new Mutation(nodeNumber, patterns[i], fromStates[i], toStates[i]);
        }
        return mutations;
    }

    /**
     * Accumulates a branch's events: the hits of each site must be added together, in path order.
     */
    public static final class Builder {
        private final boolean timed;
        private int[] subPatterns = new int[16];
        private byte[] subFrom = new byte[16];
        private byte[] subTo = new byte[16];
        private double[] subTimes;
        private int subCount = 0;
        private int[] ambPatterns = new int[16];
        private byte[] ambFrom = new byte[16];
        private int ambCount = 0;

        public Builder(boolean timed) {
            this.timed = timed;
            this.subTimes = timed ? new double[16] : null;
        }

        public void addSubstitution(int pattern, int from, int to, double time) {
            if (subCount == subPatterns.length) {
                int n = 2 * subCount;
                subPatterns = Arrays.copyOf(subPatterns, n);
                subFrom = Arrays.copyOf(subFrom, n);
                subTo = Arrays.copyOf(subTo, n);
                if (timed) {
                    subTimes = Arrays.copyOf(subTimes, n);
                }
            }
            subPatterns[subCount] = pattern;
            subFrom[subCount] = (byte) from;
            subTo[subCount] = (byte) to;
            if (timed) {
                subTimes[subCount] = time;
            }
            subCount++;
        }

        public void addAmbiguity(int pattern, int from) {
            if (ambCount == ambPatterns.length) {
                ambPatterns = Arrays.copyOf(ambPatterns, 2 * ambCount);
                ambFrom = Arrays.copyOf(ambFrom, 2 * ambCount);
            }
            ambPatterns[ambCount] = pattern;
            ambFrom[ambCount++] = (byte) from;
        }

        public BranchMutations build() {
            int size = subCount + ambCount;
            if (size == 0) {
                return EMPTY;
            }
            // events are kept in pattern order (several hits at a site stay in the order they were added)
            sortByPattern(subPatterns, subCount, subFrom, subTo, subTimes);
            sortByPattern(ambPatterns, ambCount, ambFrom, null, null);
            int[] patterns = new int[size];
            byte[] from = new byte[size];
            byte[] to = new byte[size];
            System.arraycopy(subPatterns, 0, patterns, 0, subCount);
            System.arraycopy(subFrom, 0, from, 0, subCount);
            System.arraycopy(subTo, 0, to, 0, subCount);
            System.arraycopy(ambPatterns, 0, patterns, subCount, ambCount);
            System.arraycopy(ambFrom, 0, from, subCount, ambCount);
            Arrays.fill(to, subCount, size, (byte) -1);
            double[] times = timed ? Arrays.copyOf(subTimes, size) : null;
            return new BranchMutations(patterns, from, to, subCount, times);
        }

        /**
         * Sorts the first count entries by pattern, keeping the order of equal patterns.
         */
        private static void sortByPattern(int[] patterns, int count, byte[] from, byte[] to, double[] times) {
            boolean sorted = true;
            for (int i = 1; i < count && sorted; i++) {
                sorted = patterns[i - 1] <= patterns[i];
            }
            if (sorted) {
                return;
            }
            long[] keys = new long[count];
            for (int i = 0; i < count; i++) {
                keys[i] = ((long) patterns[i] << 32) | i;
            }
            Arrays.sort(keys);
            int[] newPatterns = new int[count];
            byte[] newFrom = new byte[count];
            byte[] newTo = to == null ? null : new byte[count];
            double[] newTimes = times == null ? null : new double[count];
            for (int i = 0; i < count; i++) {
                int k = (int) keys[i];
                newPatterns[i] = patterns[k];
                newFrom[i] = from[k];
                if (to != null) {
                    newTo[i] = to[k];
                }
                if (times != null) {
                    newTimes[i] = times[k];
                }
            }
            System.arraycopy(newPatterns, 0, patterns, 0, count);
            System.arraycopy(newFrom, 0, from, 0, count);
            if (to != null) {
                System.arraycopy(newTo, 0, to, 0, count);
            }
            if (times != null) {
                System.arraycopy(newTimes, 0, times, 0, count);
            }
        }
    }
}
