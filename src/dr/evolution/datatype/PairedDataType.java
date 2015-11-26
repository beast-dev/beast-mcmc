/*
 * PairedDataType.java
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
 * @version $Id: PairedDataType.java,v 1.6 2005/06/30 09:54:51 rambaut Exp $
 *
 * @author Alexei Drummond
 */
 public class PairedDataType extends DataType {

	public static final String DESCRIPTION = "pairedDataType";
	public static final PairedDataType PAIRED_NUCLEOTIDES = new PairedDataType(Nucleotides.INSTANCE);
	public static final PairedDataType PAIRED_AMINO_ACIDS = new PairedDataType(AminoAcids.INSTANCE);
	public static final PairedDataType PAIRED_CODONS = new PairedDataType(Codons.UNIVERSAL);

	public PairedDataType(DataType baseDataType) {
		this.baseDataType = baseDataType;

		stateCount = baseDataType.getStateCount() * baseDataType.getStateCount();
		ambiguousStateCount = stateCount + 2;
	}

	public final int getState(int state1, int state2) {
        if (baseDataType.isAmbiguousState(state1) || baseDataType.isAmbiguousState(state2)) {
            return getUnknownState();
        }
		return (state1 * baseDataType.getStateCount()) + state2;
	}

	public final int getState(char c1, char c2) {
		return getState(baseDataType.getState(c1), baseDataType.getState(c2));
	}

	public final int getFirstState(int state) {
		return state / baseDataType.getStateCount();
	}

	public final int getSecondState(int state) {
		return state % baseDataType.getStateCount();
	}

    @Override
    public char[] getValidChars() {
        return null;
    }

    public final int getState(char c)
	{
		throw new IllegalArgumentException("Paired datatype cannot be expressed as char");
	}

	public int getUnknownState() {
		return stateCount;
	}

	public int getGapState() {
		return stateCount + 1;
	}

	public final char getChar(int state) {
		throw new IllegalArgumentException("Paired datatype cannot be expressed as char");
	}

	public final String getTriplet(int state) {
		throw new IllegalArgumentException("Paired datatype cannot be expressed as triplets");
	}

	public final int[] getTripletStates(int state) {
		throw new IllegalArgumentException("Paired datatype cannot be expressed as triplets");
	}

	public String getDescription() {
		return DESCRIPTION;
	}

	public int getType() {
		return -1;
	}

	// Private members

	private DataType baseDataType;

}
