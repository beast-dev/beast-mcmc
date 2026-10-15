/*
 * Mutation.java
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
 * An event on the branch above a node at a single site (pattern): either a substitution from one
 * state to another or, if flagged as an ambiguity, the point at which the site becomes ambiguous
 * (ambiguous, e.g. a gap or N, in every tip below the node). An ambiguity has no to-state, and the site
 * contributes nothing to the likelihood on that branch or in the subtree below it.
 *
 * Mutation times are not stored: the MAT likelihood integrates them out analytically
 * over the length of the branch.
 */
public final class Mutation {
    private final int nodeNumber;
    private final int pattern;
    private final int fromState;
    private final int toState;
    private final boolean ambiguity;

    public Mutation(int nodeNumber, int pattern, int fromState, int toState) {
        this(nodeNumber, pattern, fromState, toState, false);
    }

    private Mutation(int nodeNumber, int pattern, int fromState, int toState, boolean ambiguity) {
        this.nodeNumber = nodeNumber;
        this.pattern = pattern;
        this.fromState = fromState;
        this.toState = toState;
        this.ambiguity = ambiguity;
    }

    /**
     * @return an event marking that the site becomes ambiguous on the branch above the node
     */
    public static Mutation ambiguity(int nodeNumber, int pattern, int fromState) {
        return new Mutation(nodeNumber, pattern, fromState, -1, true);
    }

    /**
     * @return the number of the node at the bottom of the branch carrying this mutation
     */
    public int getNodeNumber() {
        return nodeNumber;
    }

    public int getPattern() {
        return pattern;
    }

    /**
     * @return the state at the parent node, immediately above the event
     */
    public int getFromState() {
        return fromState;
    }

    /**
     * @return the state below the event, or -1 for an ambiguity
     */
    public int getToState() {
        return toState;
    }

    public boolean isAmbiguity() {
        return ambiguity;
    }

    @Override
    public String toString() {
        return "node " + nodeNumber + ": " + fromState + (pattern + 1) + (ambiguity ? "?" : String.valueOf(toState));
    }
}
