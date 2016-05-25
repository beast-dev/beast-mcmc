/*
 * MutationDeathType.java
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

/**
 * Package: MutationDeathType
 * Description:
 * <p/>
 * <p/>
 * Created by
 * Alexander V. Alekseyenko (alexander.alekseyenko@gmail.com)
 * Date: Mar 5, 2008
 * Time: 1:09:40 PM
 */
public class MutationDeathType extends DataType {
    protected static String DESCRIPTION = "MutationDeathType";

    protected static int UNKNOWN_STATE = 128;
    protected int[][] codes;
    protected char[] stateCodes;
    public int DEATHSTATE;

    public MutationDeathType(char deathCode, char extantCode) { // Constructor for pure death type
        super();
        initialize_internals();
        codes[extantCode] = new int[]{0};
        stateCodes[0] = extantCode;

        codes[deathCode] = new int[]{1};
        stateCodes[1] = deathCode;
        stateCount = 2;
        DEATHSTATE = 1;
        ambiguousStateCount = 0;
    }

    public MutationDeathType(DataType x, char deathCode) { // constructor for extention type
        super();
        int i;
        char stateCode;

        initialize_internals();
        for (i = 0; i < x.getStateCount(); ++i) {  /* Copy unique codes */
            stateCode = x.getCode(i).charAt(0);
            this.codes[stateCode] = new int[]{i};
            stateCodes[i] = stateCode;
        }
        this.codes[deathCode] = new int[]{i}; /* Append the state space with the death state */
        stateCodes[i] = deathCode;
        DEATHSTATE = i;
        stateCount = i + 1;

        for (i = 0; i < 128; ++i) {
            int state = x.getState((char) i);
            if (state > 0 && state < 128 && x.isAmbiguousState(state) && i != deathCode) {
                if (!x.isUnknownState(state)) {
                    int[] states = x.getStates(state);
                    this.codes[i] = new int[states.length];
                    System.arraycopy(states, 0, this.codes[i], 0, states.length);
                }
            }

        }
        ambiguousStateCount = x.getAmbiguousStateCount() + 1;
    }

    public void addAmbiguity(char ambiguityCode, String s) {
        if (s.length() == 0) {
            this.codes[ambiguityCode] = new int[stateCount];
            for (int i = 0; i < stateCount; ++i) {
                this.codes[ambiguityCode][i] = i;
            }
        } else {
            this.codes[ambiguityCode] = new int[s.length()];
            for (int i = 0; i < s.length(); ++i) {
                this.codes[ambiguityCode][i] = getState(s.charAt(i));
            }
        }
        ambiguousStateCount += 1;
    }

    private void initialize_internals() {
        this.codes = new int[128][]; /* stores states (w/ ambiguities) corresponding to codes */
        this.stateCodes = new char[128]; /* Stores characters corresponding to unique state codes*/
    }

    @Override
    public char[] getValidChars() {
        return null;
    }

    /**
     * Get state corresponding to a character
     *
     * @param c character
     * @return state
     */
    public int getState(char c) {
        if (codes[c] != null && codes[c].length == 1) {
            return codes[c][0];
        } else {
            return c;
        }
    }

    /**
     * Get state corresponding to an unknown
     *
     * @return state
     */
    public int getUnknownState() {
        return UNKNOWN_STATE;
    }


    /**
     * Get character corresponding to a given state
     *
     * @param state state
     *              <p/>
     *              return corresponding character
     */
    public char getChar(int state) {
        if (state < stateCount)
            return stateCodes[state];
        return super.getChar(state);
    }

    /**
     * returns an array containing the non-ambiguous states that this state represents.
     */
    public int[] getStates(int state) {
        if (state < stateCount)
            return codes[stateCodes[state]];
        else
            return codes[state];
    }

    /**
     * returns anarray of indicatiors for non-ambiguous states that this state represents.
     */
    public boolean[] getStateSet(int state) {
        boolean[] stateSet = new boolean[stateCount];
        int states[];
        int i;
        for (i = 0; i < stateCount; ++i)
            stateSet[i] = false;
        states = getStates(state);

        for (i = 0; states != null && i < states.length; ++i) {
            stateSet[states[i]] = true;
        }
        return stateSet;
    }

    /**
     * @return true if this character is an ambiguous state
     */
    public boolean isAmbiguousChar(char c) {
        return codes[c] != null && codes[c].length > 1;
    }

    /**
     * @return true if this character is a gap
     */
    public boolean isUnknownChar(char c) {
        return codes[c] == null;
    }

    /**
     * returns true if this state is an ambiguous state.
     */
    public boolean isAmbiguousState(int state) {
        return state >= stateCount;
    }

    /**
     * @return true if this state is an unknown state
     */
    public boolean isUnknownState(int state) {
        return state >= stateCount && codes[state] == null;
    }

    public String getDescription() {
        return DESCRIPTION;
    }

    public int getType() {
        return 314;
    }
}
