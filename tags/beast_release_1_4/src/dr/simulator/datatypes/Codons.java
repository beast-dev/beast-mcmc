/*
 * Codons.java
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
 * Implements DataType for codons
 *
 * Codons have tree different representations:
 * 	State numbers - 0-63 + 64, 65 as unknown and gap
 *  State chars - the above converted into chars, starting at 'A'
 *		and '?' + '-' for unknown and gap
 *	Strings or chars of three nucleotide characters
 *
 * @version $Id: Codons.java,v 1.1 2005/04/28 22:45:25 rambaut Exp $
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 */
public class Codons extends DataType {

	/**
	 * Name of data type. For XML and human reading of data type.
	 */
	public static final String DESCRIPTION = "codon";

	public static final Codons UNIVERSAL = new Codons(GeneticCode.UNIVERSAL);
	public static final Codons VERTEBRATE_MT = new Codons(GeneticCode.VERTEBRATE_MT);
	public static final Codons YEAST = new Codons(GeneticCode.YEAST);
	public static final Codons MOLD_PROTOZOAN_MT = new Codons(GeneticCode.MOLD_PROTOZOAN_MT);
	public static final Codons MYCOPLASMA = new Codons(GeneticCode.MYCOPLASMA);
	public static final Codons INVERTEBRATE_MT = new Codons(GeneticCode.INVERTEBRATE_MT);
	public static final Codons CILIATE = new Codons(GeneticCode.CILIATE);
	public static final Codons ECHINODERM_MT = new Codons(GeneticCode.ECHINODERM_MT);
	public static final Codons EUPLOTID_NUC = new Codons(GeneticCode.EUPLOTID_NUC);
	public static final Codons BACTERIAL = new Codons(GeneticCode.BACTERIAL);
	public static final Codons ALT_YEAST = new Codons(GeneticCode.ALT_YEAST);
	public static final Codons ASCIDIAN_MT = new Codons(GeneticCode.ASCIDIAN_MT);
	public static final Codons FLATWORM_MT = new Codons(GeneticCode.FLATWORM_MT);
	public static final Codons BLEPHARISMA_NUC = new Codons(GeneticCode.BLEPHARISMA_NUC);

	public static final int STATE_COUNT = 64;
	public static final int AMBIGUOUS_STATE_COUNT = 66;

	public static final int UNKNOWN_STATE = 64;
	public static final int GAP_STATE = 65;

	public static final String[] CODON_TRIPLETS = {
		"AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG","ACT",
		"AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG","ATT",
		"CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG","CCT",
		"CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG","CTT",
		"GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG","GCT",
		"GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG","GTT",
		"TAA", "TAC", "TAG", "TAT", "TCA", "TCC", "TCG","TCT",
		"TGA", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG","TTT",
		"???", "---"
	};

	/**
	 * Private constructor - DEFAULT_INSTANCE provides the only instance
	 */
	private Codons(GeneticCode geneticCode) {
		this.geneticCode = geneticCode;
	}

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
	public final byte getState(char c)
	{
		throw new IllegalArgumentException("Codons datatype cannot be expressed as char");
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
	 * Get state corresponding to a nucleotide triplet
	 *
	 * @param nuc1 the codon triplet as chars
	 * @param nuc2 the codon triplet as chars
	 * @param nuc3 the codon triplet as chars
	 *
	 * @return state
	 */
	public final byte getState(char nuc1, char nuc2, char nuc3)
	{
		return getState(DataType.NUCLEOTIDES.getState(nuc1),
						DataType.NUCLEOTIDES.getState(nuc2),
						DataType.NUCLEOTIDES.getState(nuc3));
	}

	/**
	 * Get state corresponding to a nucleotide triplet
	 *
	 * @param nuc1 the codon triplet as states
	 * @param nuc2 the codon triplet as states
	 * @param nuc3 the codon triplet as states
	 *
	 * @return state
	 */
	public final byte getState(byte nuc1, byte nuc2, byte nuc3)
	{
		if (DataType.NUCLEOTIDES.isGapState(nuc1) ||
			DataType.NUCLEOTIDES.isGapState(nuc2) ||
			DataType.NUCLEOTIDES.isGapState(nuc3)) {
			return GAP_STATE;
		}

		if (DataType.NUCLEOTIDES.isAmbiguousState(nuc1) ||
			DataType.NUCLEOTIDES.isAmbiguousState(nuc2) ||
			DataType.NUCLEOTIDES.isAmbiguousState(nuc3)) {
			return UNKNOWN_STATE;
		}

		return (byte)((nuc1 * 16) + (nuc2 * 4) + nuc3);
	}

	/**
	 * Get character corresponding to a given state
	 *
	 * @param state state
	 *
	 * return corresponding character
	 */
	public final char getChar(byte state)
	{
		throw new IllegalArgumentException("Codons datatype cannot be expressed as char");
	}

	/**
	 * Get triplet string corresponding to a given state
	 *
	 * @param state state
	 *
	 * return corresponding triplet string
	 */
	public final String getTriplet(byte state)
	{
		return CODON_TRIPLETS[state];
	}

	public byte[] getStates(byte state) {


		if (!isAmbiguousState(state)) {
			return new byte[] { state };
		}

		byte[] states = new byte[STATE_COUNT];
		for (byte i = 0; i < STATE_COUNT; i++) {
			states[i] = i;
		}
		return states;
	}

	public boolean[] getStateSet(byte state) {

		boolean[] stateSet = new boolean[STATE_COUNT];

		if (!isAmbiguousState(state)) {
			for (int i = 0; i < STATE_COUNT; i++) {
				stateSet[i] = false;
			}
			stateSet[state] = true;
		} else {
			for (int i = 0; i < STATE_COUNT; i++) {
				stateSet[i] = true;
			}
		}

		return stateSet;
	}

	/**
	 * Get an array of three nucleotide states making this codon state
	 *
	 * @param state state
	 *
	 * return corresponding triplet string
	 */
	public final byte[] getTripletStates(byte state)
	{
		byte[] triplet = new byte[3];

		triplet[0] = DataType.NUCLEOTIDES.getState(CODON_TRIPLETS[state].charAt(0));
		triplet[1] = DataType.NUCLEOTIDES.getState(CODON_TRIPLETS[state].charAt(1));
		triplet[2] = DataType.NUCLEOTIDES.getState(CODON_TRIPLETS[state].charAt(2));

		return triplet;
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
	 * @return the genetic code
	 */
	public GeneticCode getGeneticCode() {
		return geneticCode;
	}

    /**
     * Takes non-canonical state and returns true if it represents stop codon
     * @param state
     * @return
     */
    public final boolean isStopCodon(byte state) {
        return geneticCode.isStopCodon(state);
    }

	// Private members

	private GeneticCode geneticCode;

}
