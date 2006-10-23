/*
 * HiddenNucleotides.java
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
 * implements DataType for nucleotides with hidden state(s)
 *
 * @version $Id: HiddenNucleotides.java,v 1.4 2005/05/24 20:25:56 rambaut Exp $
 *
 * @author Alexei Drummond
 */
 public class HiddenNucleotides extends DataType {
	
	public static final String DESCRIPTION = "hiddenNucleotide";
	public static final HiddenNucleotides INSTANCE = new HiddenNucleotides(2);

    public HiddenNucleotides(int numHiddenStates) {
		hiddenStateCount = numHiddenStates;

        // @todo generalize this data type to allow more than two hidden states!
        if (hiddenStateCount != 2) throw new IllegalArgumentException("Haven't implemented more than 2 yet!");
		stateCount = 4 * hiddenStateCount;
		ambiguousStateCount = stateCount + 6;
	}

	/**
	 * Get state corresponding to a character
	 *
	 * @param c character
	 *
	 * @return state
	 */
	public int getState(char c) {
	
		switch (c) {
			case '0': return 0;
			case '1': return 1;
			case '2': return 2;
			case '3': return 3;
			case '4': return 4;
			case '5': return 5;
			case '6': return 6;
			case '7': return 7;
			case 'a': case 'A': return 8;
			case 'c': case 'C': return 9;
			case 'g': case 'G': return 10;
			case 't': case 'T': case 'u': case 'U': return 11;
		}
		return getGapState();
	}

	/**
	 * Get character corresponding to a given state
	 *
	 * @param state state
	 *
	 * return corresponding character
	 */
	public char getChar(int state) {
		
		switch (state) {
			case 0: return '0';
			case 1: return '1';
			case 2: return '2';
			case 3: return '3';
			case 4: return '4';
			case 5: return '5';
			case 6: return '6';
			case 7: return '7';
			case 8: return 'A';
			case 9: return 'C';
			case 10: return 'G';
			case 11: return 'T';
		}
		
		return '-';
	}
	
	public int[] getStates(int state) {
		if (state >= 8 && state <=11) {
			int[] states = new int[2];
			states[0] = state % 4;
			states[1] = state % 4 + 4;
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
		return (state >= 8);
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

	public int getType() { return 999; }
		
	public String getDescription() { return "Hidden-state Nucleotides"; }	
		
	private int hiddenStateCount;

    static {
        registerDataType(HiddenNucleotides.DESCRIPTION, HiddenNucleotides.INSTANCE);
    }
}
