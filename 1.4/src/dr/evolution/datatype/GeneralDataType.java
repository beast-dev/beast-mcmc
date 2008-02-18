/*
 * GeneralDataType.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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

package dr.evolution.datatype;

import dr.util.Identifiable;

/**
 * Implements a general DataType for any number of states
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: GeneralDataType.java,v 1.11 2005/05/24 20:25:56 rambaut Exp $
 */
public class GeneralDataType extends DataType implements Identifiable {

    public static final String GENERAL_DATA_TYPE = "generalDataType";
    public static final String DESCRIPTION = GENERAL_DATA_TYPE;
    public static final int TYPE = GENERAL;

    /**
     * Unlike the other standard data types, this general one has a public
     * constructor so multiple instances can be created.
     *
     * @param stateCodes the codes of the states
     */
    public GeneralDataType(char[] stateCodes) {
        this.stateCount = stateCodes.length;

        for (int i = 0; i < 128; i++) stateNumbers[i] = -1;

        states = new State[stateCount];

        for (int i = 0; i < stateCount; i++) {
            stateNumbers[(int) stateCodes[i]] = i;
            states[i] = new State(i, stateCodes[i]);
        }

        this.ambiguousStateCount = 0;

    }

    /**
     * Add an alias (a state code that represents a particular state).
     * Note that all this does is put an extra entry in the stateNumbers
     * array.
     *
     * @param code  a character that represents the state
     * @param state the state number
     */
    public void addAlias(char code, int state) {
        stateNumbers[(int) code] = state;
    }

    /**
     * Add an alias (a state code that represents a particular state).
     *
     * @param code            a character that represents the state
     * @param ambiguousStates the set of states that this char code refers to.
     */
    public void addAmbiguity(char code, int[] ambiguousStates) {
        int n = ambiguousStateCount + stateCount;

        ambiguousStateCount++;
        stateNumbers[(int) code] = n;

        State[] newStates = new State[n + 1];

        System.arraycopy(states, 0, newStates, 0, n);
        newStates[n] = new State(n, code, ambiguousStates);
        states = newStates;
    }

    /**
     * Get state corresponding to a character
     *
     * @param c character
     * @return state
     */
    public int getState(char c) {

        return stateNumbers[(int) c];
    }

    /**
     * Get state corresponding to an unknown
     *
     * @return state
     */
    public int getUnknownState() {
        return stateCount + ambiguousStateCount;
    }

    /**
     * Get state corresponding to a gap
     *
     * @return state
     */
    public int getGapState() {
        return getUnknownState();
    }

    /**
     * Get character corresponding to a given state
     *
     * @param state state
     *              <p/>
     *              return corresponding character
     */
    public char getChar(int state) {
        return states[state].code;
    }

    /**
     * returns an array containing the non-ambiguous states
     * that this state represents.
     */
    public int[] getStates(int state) {

        return states[state].ambiguities;
    }

    /**
     * returns an array containing the non-ambiguous states that this state represents.
     */
    public boolean[] getStateSet(int state) {

        boolean[] stateSet = new boolean[stateCount];
        for (int i = 0; i < stateCount; i++)
            stateSet[i] = false;

        for (int i = 0, n = states[state].ambiguities.length; i < n; i++) {
            stateSet[states[state].ambiguities[i]] = true;
        }

        return stateSet;
    }

    /**
     * description of data type
     *
     * @return string describing the data type
     */
    public String getDescription() {
        if (id != null) {
            return id;
        } else {
            return DESCRIPTION;
        }
    }

    /**
     * type of data type
     *
     * @return integer code for the data type
     */
    public int getType() {
        return TYPE;
    }

    // **************************************************************
    // Identifiable IMPLEMENTATION
    // **************************************************************

    private String id = null;

    public void setId(String id) {
        this.id = id;
    }

    public String getId() {
        return id;
    }

    private int[] stateNumbers = new int[128];
    private State[] states;

    private class State {
        int number;
        char code;

        int[] ambiguities;

        State(int number, char code) {
            this.number = number;
            this.code = code;
            this.ambiguities = new int[]{number};
        }

        State(int number, char code, int[] ambiguities) {
            this.number = number;
            this.code = code;
            this.ambiguities = ambiguities;
        }
    }

}
