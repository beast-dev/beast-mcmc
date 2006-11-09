/*
 * Nucleotides.java
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

package dr.simulator.datatypes;

/**
 * implements DataType for nucleotides with ambiguous characters
 *
 * @version $Id: Nucleotides.java,v 1.1 2005/04/28 22:45:25 rambaut Exp $
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 */
public class Nucleotides extends DataType {

	/**
	 * Name of data type. For XML and human reading of data type.
	 */
	public static final String DESCRIPTION = "nucleotide";
	public static final Nucleotides INSTANCE = new Nucleotides();

	public static final int STATE_COUNT = 4;
	public static final int AMBIGUOUS_STATE_COUNT = 18;

	public static final byte A_STATE = 0;
	public static final byte C_STATE = 1;
	public static final byte G_STATE = 2;
	public static final byte UT_STATE = 3;

	public static final byte UNKNOWN_STATE = 16;
	public static final byte GAP_STATE = 17;

	/**
	 * A table to translate state numbers (0-17) into character codes
	 */
	public static final char[] NUCLEOTIDE_CHARS =
		{ 'A','C','G','T','U','K','M','R','S','W','Y','B','D','H','V','N', UNKNOWN_CHARACTER,GAP_CHARACTER};

	/**
	 * This table maps nucleotide characters into state codes (0-17)
	 * Nucleotides go ACGTURYMWSKBDHVN?-", Other letters are mapped to ?.
	 * ? and - are mapped to themselves. All other chars are mapped to -.
	 */
	public static final byte NUCLEOTIDE_STATES[] = {
		17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,	// 0-15
		17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,	// 16-31
	//                                          -
		17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,	// 32-47
	//                                                ?
		17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,16,	// 48-63
	//	    A  B  C  D  e  f  G  H  i  j  K  l  M  N  o
		17, 0,11, 1,12,16,16, 2,13,16,16,10,16, 7,15,16,	// 64-79
	//	 p  q  R  S  T  U  V  W  x  Y  z
		16,16, 5, 9, 3, 4,14, 8,16, 6,16,17,17,17,17,17,	// 80-95
	//	    A  B  C  D  e  f  G  H  i  j  K  l  M  N  o
		17, 0,11, 1,12,16,16, 2,13,16,16,10,16, 7,15,16,	// 96-111
	//	 p  q  R  S  T  U  V  W  x  Y  z
		16,16, 5, 9, 3, 4,14, 8,16, 6,16,17,17,17,17,17		// 112-127
	};

	/**
	 * A table to map state numbers (0-17) to their ambiguities
	 */
	public static final String[] NUCLEOTIDE_AMBIGUITIES = {
	//	 A    C	   G    T    U    R     Y     M     W     S     K
		"A", "C", "G", "T", "T", "AG", "CT", "AC", "AT", "CG", "GT",
	//   B      D      H      V      N       ?       -
		"CGT", "AGT", "ACT", "ACG", "ACGT", "ACGT", "ACGT"
	};

	/**
	 * Private constructor - DEFAULT_INSTANCE provides the only instance
	 */
	private Nucleotides() { }

	public int getStateCount() {
		return STATE_COUNT;
	}

	public int getAmbiguousStateCount() {
		return AMBIGUOUS_STATE_COUNT;
	}

	/**
	 * Get state corresponding to a character
	 *
	 * @param c character
	 *
	 * @return state
	 */
	public byte getState(char c) {
		return NUCLEOTIDE_STATES[c];
	}

	/**
	 * Get state corresponding to an unknown
	 *
	 * @return state
	 */
	public byte getUnknownState() {
		return UNKNOWN_STATE;
	}

	/**
	 * Get state corresponding to a gap
	 *
	 * @return state
	 */
	public byte getGapState() {
		return GAP_STATE;
	}

	/**
	 * Get character corresponding to a given state
	 *
	 * @param state state
	 *
	 * return corresponding character
	 */
	public char getChar(byte state) {
		return NUCLEOTIDE_CHARS[state];
	}

	public String getTriplet(byte state) {
		return " " + getChar(state) + " ";
	}

	/**
	 * returns an array containing the non-ambiguous states
	 * that this state represents.
	 */
	public byte[] getStates(byte state) {

		String stateString = NUCLEOTIDE_AMBIGUITIES[state];
		byte[] states = new byte[stateString.length()];
		for (int i = 0; i < stateString.length(); i++) {
			states[i] = getState(stateString.charAt(i));
		}

		return states;
	}

	/**
	 * returns an array containing the non-ambiguous states that this state represents.
	 */
	public boolean[] getStateSet(byte state) {

		boolean[] stateSet = new boolean[STATE_COUNT];
		for (int i = 0; i < STATE_COUNT; i++)
			stateSet[i] = false;

		int len = NUCLEOTIDE_AMBIGUITIES[state].length();
		for (int i = 0; i < len; i++)
			stateSet[getState(NUCLEOTIDE_AMBIGUITIES[state].charAt(i))] = true;

		return stateSet;
	}

	/**
	 * description of data type
	 *
	 * @return string describing the data type
	 */
	public String getDescription() {
		return DESCRIPTION;
	}
}
