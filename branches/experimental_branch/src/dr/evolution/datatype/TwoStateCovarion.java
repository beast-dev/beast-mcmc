/*
 * TwoStateCovarion.java
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

/**
 * @author Helen Shearman
 * @author Alexei Drummond
 */
public class TwoStateCovarion extends DataType {

    public static final String DESCRIPTION = "twoStateCovarion";
    public static final int TYPE = COVARION;
    public static final TwoStateCovarion INSTANCE = new TwoStateCovarion();

    private TwoStateCovarion() {
        stateCount = 4;
        ambiguousStateCount = 8;
    }

    @Override
    public char[] getValidChars() {
        return null;
    }

    public int getState(char c) {
        switch (c) {
            case 'a':
                return 0;
            case 'b':
                return 1;
            case 'c':
                return 2;
            case 'd':
                return 3;
            case '0':
                return 4;
            case '1':
                return 5;
            case '?':
                return getUnknownState();
            case '-':
                return getGapState();
        }
        throw new IllegalArgumentException("Character " + c + " not recognised in two-state covarion datatype!");
    }

    public char getChar(int state) {

        switch (state) {
            case 0:
                return 'a';
            case 1:
                return 'b';
            case 2:
                return 'c';
            case 3:
                return 'd';
            case 4:
                return '0';
            case 5:
                return '1';
            case 6:
                return '?';
            case 7:
                return '-';
        }
        throw new IllegalArgumentException("State " + state + " not recognised in two-state covarion datatype!");
    }

    /**
     * returns an array containing the non-ambiguous states
     * that this state represents.
     */
    public int[] getStates(int state) {

        if (state == 4 || state == 5) {
            int[] states = new int[2];
            states[0] = state % 2;
            states[1] = state % 2 + 2;
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
        } else if (state == 4 || state == 5) {
            for (int i = 0; i < stateCount; i++) {
                if ((i % 2) == (state % 2)) {
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

    public int getUnknownState() {
        return stateCount + 2;
    }

    public int getGapState() {
        return stateCount + 3;
    }

    public boolean isAmbiguousChar(char c) {
        return isAmbiguousState(getState(c));
    }

    public boolean isUnknownChar(char c) {
        return isUnknownState(getState(c));
    }

    public boolean isGapChar(char c) {
        return isGapState(getState(c));
    }

    public boolean isAmbiguousState(int state) {
        return (state >= 4);
    }

    public boolean isUnknownState(int state) {
        return (state == getUnknownState());
    }

    public boolean isGapState(int state) {
        return (state == getGapState());
    }


    public String getDescription() {
        return DESCRIPTION;
    }

    public int getType() {
        return TYPE;
    }

}
