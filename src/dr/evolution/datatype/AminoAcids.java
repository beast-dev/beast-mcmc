/*
 * AminoAcids.java
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
 * implements DataType for amino acids.
 *
 * @version $Id: AminoAcids.java,v 1.12 2005/06/22 14:48:19 beth Exp $
 *
 * @author Andrew Rambaut
 * @author Korbinian Strimmer
 * @author Alexei Drummond
 * @author Matthew Goode
 */
public class AminoAcids extends DataType
{
	/**
	 * Name of data type. For XML and human reading of data type.
	 */
	public static final String DESCRIPTION = "amino acid";

    /**
     * This integer is a unique identifier of this data type.
     */
    public static final int TYPE = AMINO_ACIDS;

    /**
     * The only instance of the AminoAcids class.
     */
    public static final AminoAcids INSTANCE = new AminoAcids();

    /**
     * This character represents the amino acid equivalent of a stop codon to cater for
     * situations arising from converting coding DNA to an amino acid sequence.
     */
    public static final char STOP_CHARACTER = '*';
	
	public static final int STOP_STATE = 23;

    /**
     * This state represents a amino acid residue of unknown type.
     */
    public static final int UNKNOWN_STATE = 24;

    /**
     * This state represents a gap in an amino acid sequence.
     */
    public static final int GAP_STATE = 25;

    /**
     * Unique integer identifier for the amino acid data type.
     */
	public static final int AMINOACIDS = 1;

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
	public static final int[] AMINOACID_STATES = {
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
	protected AminoAcids() {
        super();
        stateCount = 20;
        ambiguousStateCount = 26;
    }

    @Override
    public char[] getValidChars() {
        return AMINOACID_CHARS;
    }

    /**
	 * Get state corresponding to a character.
	 *
	 * @param c character
	 *
	 * @return state
	 */
	public int getState(final char c) {

		return AMINOACID_STATES[c];
	}

	/**
	 * Get state corresponding to a stop.
	 *
	 * @return state
	 */
	public static int getStopState() {

		return STOP_STATE;
	}

	/**
	 * Get state corresponding to an unknown.
	 *
	 * @return state
	 */
	public int getUnknownState() {

		return AminoAcids.UNKNOWN_STATE;
	}

	/**
	 * Get state corresponding to a gap.
	 *
	 * @return state
	 */
	public int getGapState() {
		return AminoAcids.GAP_STATE;
	}

	/**
	 * Get character corresponding to a given state.
	 *
	 * @param state state
	 *
	 * @return corresponding character
	 */
	public char getChar(final int state) {
		return AminoAcids.AMINOACID_CHARS[state];
	}

	/**
	 * Get triplet string corresponding to a given state.
	 *
	 * @param state state
	 *
	 * @return corresponding triplet string
	 */
	public String getTriplet(final int state) {

		return AminoAcids.AMINOACID_TRIPLETS[state];
	}

	/**
     * @param state the state to return the state set of. If this state is an ambiguity
     * code then the array returned will be size greater than 1.
	 * @return an array containing the non-ambiguous states
	 * that this state represents.
	 */
	public int[] getStates(final int state) {

		final String stateString = AminoAcids.AMINOACID_AMBIGUITIES[state];
		final int[] states = new int[stateString.length()];
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
	public boolean[] getStateSet(final int state) {
	
		final boolean[] stateSet = new boolean[stateCount];
		for (int i = 0; i < stateCount; i++) {
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
	 * type of data type.
	 *
	 * @return integer code for the data type
	 */
	public int getType() {
		return TYPE;
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
	public boolean isStopState(final int state) {
		return state == getStopState();
	}
}
