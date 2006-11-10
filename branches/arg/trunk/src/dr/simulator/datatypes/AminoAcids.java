/*
 * AminoAcids.java
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
 * implements DataType for amino acids.
 *
 * @version $Id: AminoAcids.java,v 1.1 2005/04/28 22:45:25 rambaut Exp $
 *
 * @author Andrew Rambaut
 * @author Korbinian Strimmer
 * @author Alexei Drummond
 * @author Matthew Goode
 */
public final class AminoAcids extends DataType
{
	/**
	 * Name of data type. For XML and human reading of data type.
	 */
	public static final String DESCRIPTION = "amino acid";

    /**
     * The only instance of the AminoAcids class.
     */
    public static final AminoAcids INSTANCE = new AminoAcids();

	public static final int STATE_COUNT = 20;
	public static final int AMBIGUOUS_STATE_COUNT = 25;

	/**
     * This character represents the amino acid equivalent of a stop codon to cater for
     * situations arising from converting coding DNA to an amino acid sequences.
     */
    public static final char STOP_CHARACTER = '*';

	public static final byte STOP_STATE = 23;

    /**
     * This state represents a amino acid residue of unknown type.
     */
    public static final byte UNKNOWN_STATE = 24;

    /**
     * This state represents a gap in an amino acid sequences.
     */
    public static final byte GAP_STATE = 25;

    /**
     * Unique integer identifier for the amino acid data type.
     */
	public static final byte AMINOACIDS = 1;

	/**
	 * A table to translate state numbers (0-25) into one letter codes.
	 */
	public static final char[] AMINOACID_CHARS= {
		'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R',
		'S','T','V','W','Y','B','Z','X',AminoAcids.STOP_CHARACTER,DataType.UNKNOWN_CHARACTER,DataType.GAP_CHARACTER
	};

	/**
	 * A table to map state numbers (0-25) to their three letter codes.
	 */
	private static final String[] AMINOACID_TRIPLETS = {
	//		A		C		D		E		F		G		H		I		K
		  "Ala",  "Cys",  "Asp",  "Glu",  "Phe",  "Gly",  "His",  "Ile",  "Lys",
	//		L		M		N		P		Q		R		S		T		V
		  "Leu",  "Met",  "Asn",  "Pro",  "Gln",  "Arg",  "Ser",  "Thr",  "Val",
	//		W		Y		B		Z		X		*		?		-
		  "Trp",  "Tyr",  "Asx",  "Glx",  " X ",  " * ",  " ? ",  " - "
	};

	/**
	 * This table maps amino acid characters into state codes (0-25).
	 * Amino Acids go ACDEFGHIKLMNPQRSTVWYBZX*?-,
	 * Other letters; j, o, and u are mapped to ?
	 * *, ? and - are mapped to themselves
	 * All other chars are mapped to -
	 */
	public static final byte[] AMINOACID_STATES = {
		25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,	// 0-15
		25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,	// 16-31
	//                                 *        -
		25,25,25,25,25,25,25,25,25,25,23,25,25,25,25,25,	// 32-47
	//                                                ?
		25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,24,	// 48-63
	//		A  B  C  D  E  F  G  H  I  j  K  L  M  N  o
		25, 0,20, 1, 2, 3, 4, 5, 6, 7,24, 8, 9,10,11,24,	// 64-79
	//	 P  Q  R  S  T  u  V  W  X  Y  Z
		12,13,14,15,16,24,17,18,22,19,21,25,25,25,25,25,	// 80-95
	//		A  B  C  D  E  F  G  H  I  j  K  L  M  N  o
		25, 0,20, 1, 2, 3, 4, 5, 6, 7,24, 8, 9,10,11,24,	// 96-111
	//	 P  Q  R  S  T  u  V  W  X  Y  Z
		12,13,14,15,16,24,17,18,22,19,21,25,25,25,25,25		// 112-127
	};

	/**
	 * A table to map state numbers (0-25) to their ambiguities.
	 */
	private static final String[] AMINOACID_AMBIGUITIES = {
	//	   A	C	 D	  E	   F	G	 H	  I	   K
		  "A", "C", "D", "E", "F", "G", "H", "I", "K",
	//	   L	M	 N	  P	   Q	R	 S	  T	   V
		  "L", "M", "N", "P", "Q", "R", "S", "T", "V",
	//	   W	Y	 B	   Z
		  "W", "Y", "DN", "EQ",
	//	   X					   *	?						-
		  "ACDEFGHIKLMNPQRSTVWY", "*", "ACDEFGHIKLMNPQRSTVWY", "ACDEFGHIKLMNPQRSTVWY"
	};

	/**
	 * Private constructor - DEFAULT_INSTANCE provides the only instance.
	 */
	private AminoAcids() {
        super();
    }

	public int getStateCount() {
		return STATE_COUNT;
	}

	public int getAmbiguousStateCount() {
		return AMBIGUOUS_STATE_COUNT;
	}

	/**
	 * Get state corresponding to a character.
	 *
	 * @param c character
	 *
	 * @return state
	 */
	public byte getState(final char c) {

		return AMINOACID_STATES[c];
	}

	/**
	 * Get state corresponding to a stop.
	 *
	 * @return state
	 */
	public static byte getStopState() {

		return STOP_STATE;
	}

	/**
	 * Get state corresponding to an unknown.
	 *
	 * @return state
	 */
	public byte getUnknownState() {

		return AminoAcids.UNKNOWN_STATE;
	}

	/**
	 * Get state corresponding to a gap.
	 *
	 * @return state
	 */
	public byte getGapState() {
		return AminoAcids.GAP_STATE;
	}

	/**
	 * Get character corresponding to a given state.
	 *
	 * @param state state
	 *
	 * @return corresponding character
	 */
	public char getChar(final byte state) {
		return AminoAcids.AMINOACID_CHARS[state];
	}

	/**
	 * Get triplet string corresponding to a given state.
	 *
	 * @param state state
	 *
	 * @return corresponding triplet string
	 */
	public String getTriplet(final byte state) {

		return AminoAcids.AMINOACID_TRIPLETS[state];
	}

	/**
     * @param state the state to return the state set of. If this state is an ambiguity
     * code then the array returned will be size greater than 1.
	 * @return an array containing the non-ambiguous states
	 * that this state represents.
	 */
	public byte[] getStates(final byte state) {

		final String stateString = AminoAcids.AMINOACID_AMBIGUITIES[state];
		final byte[] states = new byte[stateString.length()];
		for (int i = 0; i < stateString.length(); i++) {
			states[i] = getState(stateString.charAt(i));
		}

		return states;
	}

	/**
	 * @param state the state to return the state set of. If this state is an ambiguity
     * code then the state set returned will have more than one true element in it.
     * @return an array containing the non-ambiguous states that this state represents.
	 */
	public boolean[] getStateSet(final byte state) {

		final boolean[] stateSet = new boolean[STATE_COUNT];
		for (int i = 0; i < STATE_COUNT; i++) {
			stateSet[i] = false;
        }

		final int len = AminoAcids.AMINOACID_AMBIGUITIES[state].length();
		for (int i = 0; i < len; i++) {
			stateSet[getState(AMINOACID_AMBIGUITIES[state].charAt(i))] = true;
        }

		return stateSet;
	}

	/**
	 * description of data type.
	 *
	 * @return string describing the data type
	 */
	public String getDescription() {
		return DESCRIPTION;
	}

	/**
	 * @param c the char that is being tested to see if its a stop char.
     * @return true if this character is a stop
	 */
	public boolean isStopChar(final char c) {
		return isStopState(getState(c));
	}

	/**
     * @param state the state that is being tested to see if its a stop state.
	 * @return true if this state is a stop.
	 */
	public boolean isStopState(final byte state) {
		return state == getStopState();
	}
}
