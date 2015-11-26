/*
 * Nucleotides.java
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
 * implements DataType for nucleotides with ambiguous characters
 *
 * @version $Id: Nucleotides.java,v 1.10 2006/08/31 14:57:24 rambaut Exp $
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 */
public class Nucleotides extends DataType {

    public static final String JC = "JC";
	public static final String F84 = "F84";
	public static final String HKY = "HKY";
	public static final String GTR = "GTR";
	
	/**
	 * Name of data type. For XML and human reading of data type.
	 */
	public static final String DESCRIPTION = "nucleotide";
	public static final int TYPE = NUCLEOTIDES;
	public static final Nucleotides INSTANCE = new Nucleotides();

    public static final int A_STATE = 0;
	public static final int C_STATE = 1;
	public static final int G_STATE = 2;
	public static final int UT_STATE = 3;

    public static final int R_STATE = 5; // A or G
    public static final int Y_STATE = 6; // C or T

	public static final int UNKNOWN_STATE = 16;
	public static final int GAP_STATE = 17;
	
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
	public static final int NUCLEOTIDE_STATES[] = {
		17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,	// 0-15
		17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,	// 16-31
	//                                          -
		17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,	// 32-47
	//                                                ?
		17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,16,	// 48-63
	//	    A  B  C  D  e  f  G  H  i  j  K  l  M  N  o
		17, 0,11, 1,12,16,16, 2,13,16,16,10,16, 7,15,16,	// 64-79
	//	 p  q  R  S  T  U  V  W  x  Y  z
		16,16, 5, 9, 3, 3,14, 8,16, 6,16,17,17,17,17,17,	// 80-95
	//	    A  B  C  D  e  f  G  H  i  j  K  l  M  N  o
		17, 0,11, 1,12,16,16, 2,13,16,16,10,16, 7,15,16,	// 96-111
	//	 p  q  R  S  T  U  V  W  x  Y  z
		16,16, 5, 9, 3, 3,14, 8,16, 6,16,17,17,17,17,17		// 112-127
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
	protected Nucleotides() {
		stateCount = 4;
		ambiguousStateCount = 18;
	}

    @Override
    public char[] getValidChars() {
        return NUCLEOTIDE_CHARS;
    }

    /**
	 * Get state corresponding to a character
	 *
	 * @param c character
	 *
	 * @return state
	 */
	public int getState(char c) {
		return NUCLEOTIDE_STATES[c];
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
	 * Get state corresponding to a gap
	 *
	 * @return state
	 */
	public int getGapState() {
		return GAP_STATE;
	}

	/**
	 * Get character corresponding to a given state
	 *
	 * @param state state
	 *
	 * return corresponding character
	 */
	public char getChar(int state) {
		return NUCLEOTIDE_CHARS[state];
	}

	/**
	 * returns an array containing the non-ambiguous states
	 * that this state represents.
	 */
	public int[] getStates(int state) {

		String stateString = NUCLEOTIDE_AMBIGUITIES[state];
		int[] states = new int[stateString.length()];
		for (int i = 0; i < stateString.length(); i++) {
			states[i] = getState(stateString.charAt(i));
		}

		return states;
	}
	
	/**
	 * returns an array containing the non-ambiguous states that this state represents.
	 */
	public boolean[] getStateSet(int state) {
	
		boolean[] stateSet = new boolean[stateCount];
		for (int i = 0; i < stateCount; i++)
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

	/**
	 * type of data type
	 *
	 * @return integer code for the data type
	 */
	public int getType() {
		return TYPE;
	}

}
