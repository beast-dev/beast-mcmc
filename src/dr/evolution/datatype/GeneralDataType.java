/*
 * GeneralDataType.java
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

package dr.evolution.datatype;

import dr.util.Identifiable;

import java.util.*;

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
    public static final GeneralDataType INSTANCE = new GeneralDataType();

    // for BEAUti trait PartitionSubstitutionModel
    public GeneralDataType() {}
    /**
     * Unlike the other standard data types, this general one has a public
     * constructor so multiple instances can be created.
     *
     * @param stateCodes the codes of the states
     */
    public GeneralDataType(final String[] stateCodes) {
        for (int i = 0; i < stateCodes.length; i++) {
            State state = new State(i, stateCodes[i]);
            states.add(state);
            stateMap.put(stateCodes[i], state);
        }
        stateCount = states.size();

        this.ambiguousStateCount = 0;

    }

    /**
     * Unlike the other standard data types, this general one has a public
     * constructor so multiple instances can be created.
     *
     * @param stateCodes the codes of the states
     */
    public GeneralDataType(final Collection<String> stateCodes) {
        int i = 0;
        for (String code : stateCodes) {
            State state = new State(i, code);
            states.add(state);
            stateMap.put(code, state);
            i++;
        }
        stateCount = states.size();

        this.ambiguousStateCount = 0;
    }

    /**
     * Add an alias (a state code that represents a particular state).
     * Note that all this does is put an extra entry in the stateNumbers
     * array.
     *
     * @param alias a string that represents the state
     * @param code the state number
     */
    public void addAlias(String alias, String code) {
        State state =stateMap.get(code);
        if (state == null) {
            throw new IllegalArgumentException("DataType doesn't contain the state, " + code);
        }
        stateMap.put(alias, state);
    }

    /**
     * Add an ambiguity (a state code that represents multiple states).
     *
     * @param code            a string that represents the state
     * @param ambiguousStates the set of states that this code refers to.
     */
    public void addAmbiguity(String code, String[] ambiguousStates) {

        int n = ambiguousStateCount + stateCount;

        int[] indices = new int[ambiguousStates.length];
        int i = 0;
        for (String stateCode : ambiguousStates) {
            State state =stateMap.get(stateCode);
            if (state == null) {
                throw new IllegalArgumentException("DataType doesn't contain the state, " + stateCode);
            }
            indices[i] = state.number;
            i++;
        }
        State state = new State(n, code, indices);
        states.add(state);
        ambiguousStateCount++;

        stateMap.put(code, state);
    }

    @Override
    public char[] getValidChars() {
        return null;
    }

    /**
     * Get state corresponding to a code
     *
     * @param code string code
     * @return state
     */
    public int getState(String code) {
        if (code.equals("?")) {
            return getUnknownState();
        }
        if (!stateMap.containsKey(code)) {
            return -1;
        }
        return stateMap.get(code).number;
    }

    /**
     * Override this function to cast to string codes...
     * @param c character
     *
     * @return
     */
    public int getState(char c) {
        return getState(String.valueOf(c));
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
     * @return corresponding code
     */
    public String getCode(int state) {
        return states.get(state).code;
    }

    /**
     * returns an array containing the non-ambiguous states
     * that this state represents.
     */
    public int[] getStates(int state) {

        return states.get(state).ambiguities;
    }

    /**
     * returns an array containing the non-ambiguous states that this state represents.
     */
    public boolean[] getStateSet(int state) {

        boolean[] stateSet = new boolean[stateCount];

        if (state < states.size()) {
            State s = states.get(state);

            for (int i = 0; i < stateCount; i++) {
                stateSet[i] = false;
            }
            for (int i = 0, n = s.ambiguities.length; i < n; i++) {
                stateSet[s.ambiguities[i]] = true;
            }
        } else if (state == states.size()) {
            for (int i = 0; i < stateCount; i++) {
                stateSet[i] = true;
            }
        } else {
            throw new IllegalArgumentException("invalid state index");
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

    private List<State> states = new ArrayList<State>();
    private Map<String, State> stateMap = new TreeMap<String, State>();

    private class State {
        int number;
        String code;

        int[] ambiguities;

        State(int number, String code) {
            this.number = number;
            this.code = code;
            this.ambiguities = new int[]{number};
        }

        State(int number, String code, int[] ambiguities) {
            this.number = number;
			this.code = code;
			this.ambiguities = ambiguities;
		}
	}

}
