/*
 * CodonTable.java
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
 * Describes a device for translating Nucleotide triplets
 * or codon indices into amino acid codes.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 *
 * @version $Id: CodonTable.java,v 1.9 2005/05/24 20:25:56 rambaut Exp $
 */

public interface CodonTable {

	/**
	 * Returns the char associated with AminoAcid represented by codonState.
	 * Note that the char is as defined by AminoAcids.java
	 * @see AminoAcids
	 * @see Codons
	 * @return state for '?' if codon unknown
	 */
	char getAminoAcidChar(int codonState);

	/**
	 * Returns the state associated with AminoAcid represented by codonState.
	 * Note that the state is as defined by AminoAcids.java
	 * @see AminoAcids
	 * @see Codons
	 * @return '?' if codon unknown
	 */
	int getAminoAcidState(int codonState);

	/**
	 * @return all the possible codons for a given amino acid
	 */
	char[][] getCodonsFromAminoAcidState(int aminoAcidState);

	/*
	 * @return all the possible codons for a given amino acid
	 */
	char[][] getCodonsFromAminoAcidChar(char aminoAcidChar);

	/*
	 * @returns three IUPAC states representing the given amino acid
	 * @note The returned array should not be altered, and implementations
	 *       should attempt to implement this as efficiently as possible
	 */
	int[] getAmbiguousCodonFromAminoAcidState(int aminoAcid);

	/**
	 * @return the codon states of terminator amino acids.
	 */
	int[] getStopCodonIndices();

	/**
	 * Returns the number of terminator amino acids.
	 */
	int getStopCodonCount();
	
	/**
	 * A simple concrete class for the standard genetic codes
	 */

}
