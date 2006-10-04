/*
 * Binary.java
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
 * implements DataType for two-state data
 *
 * @version $Id: Binary.java,v 1.1 2005/04/28 22:45:25 rambaut Exp $
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 */
public class Binary extends DataType {

    /**
	 * Name of data type. For XML and human reading of data type.
	 */
	public static final String DESCRIPTION = "binary";
	public static final Binary INSTANCE = new Binary();

	public static final int STATE_COUNT = 2;
	public static final int AMBIGUOUS_STATE_COUNT = 4;

	public static final byte ZERO_STATE = 0;
	public static final byte ONE_STATE = 1;

	public static final byte UNKNOWN_STATE = 2;
	public static final byte GAP_STATE = 3;

	/**
	 * A table to translate state numbers (0-3) into character codes
	 */
	public static final char[] BINARY_CHARS =
		{ '0','1', UNKNOWN_CHARACTER,GAP_CHARACTER};

	/**
	 * A table to map state numbers (0-3) to their ambiguities
	 */
	public static final String[] BINARY_AMBIGUITIES = {
	//	 0    1	   ?     -
		"0", "1", "01", "01"
	};

	/**
	 * Private constructor - DEFAULT_INSTANCE provides the only instance
	 */
	private Binary() { }

	public int getStateCount() {
		return STATE_COUNT;
	}

	public int getAmbiguousStateCount() {
		return AMBIGUOUS_STATE_COUNT;
	}

	// Get state corresponding to character c
	public byte getState(char c) {
		switch (c) {
 			case '0':
				return 0;
			case '1':
				return 1;
			case UNKNOWN_CHARACTER:
				return 2;
			case GAP_CHARACTER:
				return 3;
			default:
				return 2;
		}
	}

	public byte getUnknownState() {
		return UNKNOWN_STATE;
	}

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
		return BINARY_CHARS[state];
	}

	public String getTriplet(byte state) {
		return " " + getChar(state) + " ";
	}

	/**
	 * returns an array containing the non-ambiguous states
	 * that this state represents.
	 */
	public byte[] getStates(byte state) {

		String stateString = BINARY_AMBIGUITIES[state];
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
		if (state < 2) {
			stateSet[1-state] = false;
			stateSet[state] = true;
		} else {
			stateSet[0] = true;
			stateSet[1] = true;
		}

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
