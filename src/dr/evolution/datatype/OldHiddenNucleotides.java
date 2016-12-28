/*
 * HiddenNucleotides.java
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
 * implements DataType for nucleotides with hidden state(s)
 *
 * @author Alexei Drummond
 * @version $Id: HiddenNucleotides.java,v 1.4 2005/05/24 20:25:56 rambaut Exp $
 */
@Deprecated
public class OldHiddenNucleotides extends DataType {

    public static final String DESCRIPTION = "hiddenNucleotide";
    public static final OldHiddenNucleotides INSTANCE = new OldHiddenNucleotides(2);

    public OldHiddenNucleotides(int numHiddenStates) {
        hiddenClassCount = numHiddenStates;

        stateCount = 4 * hiddenClassCount;
        ambiguousStateCount = stateCount + 6;
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

        switch (c) {
            case'A':
            case'a':
                return stateCount;
            case'C':
            case'c':
                return stateCount + 1;
            case'G':
            case'g':
                return stateCount + 2;
            case'T':
            case't':
            case'U':
            case'u':
                return stateCount + 3;
            case'-':
            case'?':
                return getGapState();
            default: {
                int state = (int) c - '0';
                if (c > '?') state -= 1;
                if (c > 'A') state -= 1;
                if (c > 'C') state -= 1;
                if (c > 'G') state -= 1;
                if (c > 'T') state -= 1;
                if (c > 'U') state -= 1;
                if (c > 'a') state -= 1;
                if (c > 'g') state -= 1;
                if (c > 'c') state -= 1;
                if (c > 't') state -= 1;
                if (c > 'u') state -= 1;
                return state;
            }
        }

    }

    /**
     * Get character corresponding to a given state
     *
     * @param state state
     *              <p/>
     *              return corresponding character
     */
    public char getChar(int state) {

        if (state >= stateCount) {
            switch (state - stateCount) {
                case 0:
                    return 'A';
                case 1:
                    return 'C';
                case 2:
                    return 'G';
                case 3:
                    return 'T';
                default:
                    return '-';
            }
        } else {
            char c = (char) (state + '0');
            if (c >= '?') c += 1;
            if (c >= 'A') c += 1;
            if (c >= 'C') c += 1;
            if (c >= 'G') c += 1;
            if (c >= 'T') c += 1;
            if (c >= 'U') c += 1;
            if (c >= 'a') c += 1;
            if (c >= 'g') c += 1;
            if (c >= 'c') c += 1;
            if (c >= 't') c += 1;
            if (c >= 'u') c += 1;
            return c;
        }
    }

    public int[] getStates(int state) {

        if (state >= stateCount && state <= stateCount + 3) {
            int[] states = new int[hiddenClassCount];
            for (int i = 0; i < hiddenClassCount; i++) {
                states[i] = state % 4 + (i * 4);
            }
            return states;
        } else throw new IllegalArgumentException();
    }

    /**
     * returns an array containing the non-ambiguous states that this state represents.
     */
    public boolean[] getStateSet(int state) {

        boolean[] stateSet = new boolean[stateCount];
        for (int i = 0; i < stateCount; i++) {
            stateSet[i] = false;
        }
        if (!isAmbiguousState(state)) {
            stateSet[state] = true;
        } else if (state < (stateCount + 4)) {
            for (int i = 0; i < stateCount; i++) {
                if ((i % 4) == (state % 4)) {
                    stateSet[i] = true;
                }

            }
        } else {
            for (int i = 0; i < stateCount; i++) {
                stateSet[i] = true;
            }
        }

        return stateSet;
    }

    /**
     * Get state corresponding to an unknown
     *
     * @return state
     */
    public int getUnknownState() {
        return stateCount + 4;
    }

    /**
     * Get state corresponding to a gap
     *
     * @return state
     */
    public int getGapState() {
        return stateCount + 5;
    }

    /**
     * @return true if this character is an ambiguous state
     */
    public boolean isAmbiguousChar(char c) {
        return isAmbiguousState(getState(c));
    }

    /**
     * @return true if this character is a gap
     */
    public boolean isUnknownChar(char c) {
        return isUnknownState(getState(c));
    }

    /**
     * @return true if this character is a gap
     */
    public boolean isGapChar(char c) {
        return isGapState(getState(c));
    }

    /**
     * returns true if this state is an ambiguous state.
     */
    public boolean isAmbiguousState(int state) {
        return (state >= stateCount);
    }

    /**
     * @return true if this state is an unknown state
     */
    public boolean isUnknownState(int state) {
        return (state == getUnknownState());
    }

    /**
     * @return true if this state is a gap
     */
    public boolean isGapState(int state) {
        return (state == getGapState());
    }

    public int getType() {
        return 999;
    }

    public String getDescription() {
        return "Hidden-state Nucleotides";
    }

    private int hiddenClassCount;

    public int getHiddenClassCount() {
        return hiddenClassCount;
    }
}
