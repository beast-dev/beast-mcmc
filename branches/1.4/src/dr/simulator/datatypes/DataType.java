/*
 * DataType.java
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
 * Base class for sequences data types.
 *
 * @version $Id: DataType.java,v 1.1 2005/04/28 22:45:25 rambaut Exp $
 *
 * @author Andrew Rambaut
 */
public abstract class DataType
{
	public static final DataType NUCLEOTIDES = Nucleotides.INSTANCE;
	public static final DataType AMINO_ACIDS = AminoAcids.INSTANCE;
	public static final DataType CODONS = Codons.UNIVERSAL;
	public static final DataType BINARY = Binary.INSTANCE;

	public static final char UNKNOWN_CHARACTER = '?';
	public static final char GAP_CHARACTER = '-';

	/**
	 * guess data type suitable for a given sequences
	 *
	 * @param sequence a string of symbols representing a molecular sequences of unknown data type.
	 *
	 * @return suitable DataType object
	 */
	public static DataType guessDataType(String sequence)
	{
		// count A, C, G, T, U, N
		long numNucs = 0;
		long numChars = 0;
		long numBins = 0;
		for (int i = 0; i < sequence.length(); i++)
		{
			char c = sequence.charAt(i);
			byte s = Nucleotides.INSTANCE.getState(c);

			if (s != Nucleotides.UNKNOWN_STATE && s != Nucleotides.GAP_STATE) {
				numNucs++;
			}

			if (c != '-' && c != '?') numChars++;

			if (c == '0' || c == '1') numBins++;
		}

		if (numChars == 0) { numChars = 1; }

		// more than 85 % frequency advocates nucleotide data
		if ((double) numNucs / (double) numChars > 0.85) {
			return Nucleotides.INSTANCE;
		} else if ((double) numBins / (double) numChars > 0.2) {
			return Binary.INSTANCE;
		} else {
			return AminoAcids.INSTANCE;
		}
	}

	/**
	 * Get number of unique states
	 *
	 * @return number of unique states
	 */
	public abstract int getStateCount();

	/**
	 * Get number of states including ambiguous states
	 *
	 * @return number of ambiguous states
	 */
	public abstract int getAmbiguousStateCount();

	/**
	 * Get state corresponding to a character
	 *
	 * @param c character
	 *
	 * @return state
	 */
	public abstract byte getState(char c);

	/**
	 * Get state corresponding to an unknown
	 *
	 * @return state
	 */
	public abstract byte getUnknownState();

	/**
	 * Get state corresponding to a gap
	 *
	 * @return state
	 */
	public abstract byte getGapState();

	/**
	 * Get character corresponding to a given state
	 *
	 * @param state state
	 *
	 * return corresponding character
	 */
	public abstract char getChar(byte state);

	/**
	 * Get triplet string corresponding to a given state
	 *
	 * @param state state
	 *
	 * return corresponding triplet string
	 */
	public abstract String getTriplet(byte state);

	/**
	 * returns an array containing the non-ambiguous states that this state represents.
	 */
	public abstract byte[] getStates(byte state);

	/**
	 * returns an array containing the non-ambiguous states that this state represents.
	 */
	public abstract boolean[] getStateSet(byte state);

	/**
	 * returns the uncorrected distance between two states
	 */
	public double getObservedDistance(byte state1, byte state2)
	{
		if (!isAmbiguousState(state1) && !isAmbiguousState(state2) && state1 != state2) {
			return 1.0;
		}

		return 0.0;
	}

	/**
	 * returns the uncorrected distance between two states with full
	 * treatment of ambiguity.
	 */
	public double getObservedDistanceWithAmbiguity(byte state1, byte state2)
	{
		boolean[] stateSet1 = getStateSet(state1);
		boolean[] stateSet2 = getStateSet(state2);

		double sumMatch = 0.0;
		double sum1 = 0.0;
		double sum2 = 0.0;
		for (int i = 0; i < getStateCount(); i++) {
			if (stateSet1[i]) {
				sum1 += 1.0;
				if (stateSet1[i] == stateSet2[i]) {
					sumMatch += 1.0;
				}
			}
			if (stateSet2[i]) {
				sum2 += 1.0;
			}
		}

		double distance = (1.0 - (sumMatch / (sum1 * sum2)));

		return distance;
	}

	public String toString() {
		return getDescription();
	}

	/**
	 * description of data type
	 *
	 * @return string describing the data type
	 */
	public abstract String getDescription();

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
	public boolean isAmbiguousState(byte state) {
		return (state >= getStateCount());
	}

	/**
	 * @return true if this state is an unknown state
	 */
	public boolean isUnknownState(byte state) {
		return (state == getUnknownState());
	}

	/**
	 * @return true if this state is a gap
	 */
	public boolean isGapState(byte state) {
		return (state == getGapState());
	}

}
