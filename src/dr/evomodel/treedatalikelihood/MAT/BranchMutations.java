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

/**
 * The events on one branch, held in primitive arrays: the substitutions first, then the
 * ambiguities (sites becoming ambiguous, with a to-state of -1). Instances are immutable, so a
 * branch's events can be cached against the identity of its BranchMutations object.
 */
public final class BranchMutations {

    public static final BranchMutations EMPTY = new BranchMutations(new int[0], new byte[0], new byte[0], 0);

    private final int[] patterns;
    private final byte[] fromStates;
    private final byte[] toStates;
    private final int substitutionCount;

    BranchMutations(int[] patterns, byte[] fromStates, byte[] toStates, int substitutionCount) {
        this.patterns = patterns;
        this.fromStates = fromStates;
        this.toStates = toStates;
        this.substitutionCount = substitutionCount;
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
}
