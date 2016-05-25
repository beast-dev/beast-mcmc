/*
 * SimpleCodons.java
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
 * Implements DataType for codons
 * 
 * Codons have tree different representations: 
 * 	State numbers - 0-63 + 64, 65 as unknown and gap
 *  State chars - the above converted into chars, starting at 'A'
 *		and '?' + '-' for unknown and gap
 *	Strings or chars of three nucleotide characters
 * 
 * @version $Id: SimpleCodons.java,v 1.1 2005/06/28 11:51:50 rambaut Exp $
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 */
public class SimpleCodons extends DataType {
    
	/**
	 * Name of data type. For XML and human reading of data type.
	 */
	public static final String DESCRIPTION = "simpleCodon";
	public static final int TYPE = CODONS;
	
	public static final SimpleCodons UNIVERSAL = new SimpleCodons(GeneticCode.UNIVERSAL);
	public static final SimpleCodons VERTEBRATE_MT = new SimpleCodons(GeneticCode.VERTEBRATE_MT);
	public static final SimpleCodons YEAST = new SimpleCodons(GeneticCode.YEAST);
	public static final SimpleCodons MOLD_PROTOZOAN_MT = new SimpleCodons(GeneticCode.MOLD_PROTOZOAN_MT);
	public static final SimpleCodons MYCOPLASMA = new SimpleCodons(GeneticCode.MYCOPLASMA);
	public static final SimpleCodons INVERTEBRATE_MT = new SimpleCodons(GeneticCode.INVERTEBRATE_MT);
	public static final SimpleCodons CILIATE = new SimpleCodons(GeneticCode.CILIATE);
	public static final SimpleCodons ECHINODERM_MT = new SimpleCodons(GeneticCode.ECHINODERM_MT);
	public static final SimpleCodons EUPLOTID_NUC = new SimpleCodons(GeneticCode.EUPLOTID_NUC);
	public static final SimpleCodons BACTERIAL = new SimpleCodons(GeneticCode.BACTERIAL);
	public static final SimpleCodons ALT_YEAST = new SimpleCodons(GeneticCode.ALT_YEAST);
	public static final SimpleCodons ASCIDIAN_MT = new SimpleCodons(GeneticCode.ASCIDIAN_MT);
	public static final SimpleCodons FLATWORM_MT = new SimpleCodons(GeneticCode.FLATWORM_MT);
	public static final SimpleCodons BLEPHARISMA_NUC = new SimpleCodons(GeneticCode.BLEPHARISMA_NUC);

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
	private SimpleCodons(GeneticCode geneticCode) {
		this.geneticCode = geneticCode;
		
		stateCount = 64;
		ambiguousStateCount = 66;
	}

    @Override
    public char[] getValidChars() {
        return null;
    }

    /**
	 * Get state corresponding to a character
	 *
	 * @param c character
	 *
	 * @return state
	 */
	public final int getState(char c)
	{
		throw new IllegalArgumentException("Codons datatype cannot be expressed as char");
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
	 * Get state corresponding to a nucleotide triplet
	 *
	 * @param nuc1 the codon triplet as chars
	 * @param nuc2 the codon triplet as chars
	 * @param nuc3 the codon triplet as chars
	 *
	 * @return state
	 */
	public final int getState(char nuc1, char nuc2, char nuc3)
	{
		return getState(Nucleotides.INSTANCE.getState(nuc1), 
						Nucleotides.INSTANCE.getState(nuc2), 
						Nucleotides.INSTANCE.getState(nuc3));
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
	public final int getState(int nuc1, int nuc2, int nuc3)
	{
		if (Nucleotides.INSTANCE.isGapState(nuc1) ||
			Nucleotides.INSTANCE.isGapState(nuc2) ||
			Nucleotides.INSTANCE.isGapState(nuc3)) {
			return GAP_STATE;
		}
	
		if (Nucleotides.INSTANCE.isAmbiguousState(nuc1) ||
			Nucleotides.INSTANCE.isAmbiguousState(nuc2) ||
			Nucleotides.INSTANCE.isAmbiguousState(nuc3)) {
			return UNKNOWN_STATE;
		}

        return (nuc1 * 16) + (nuc2 * 4) + nuc3;
	}

	/**
	 * Get character corresponding to a given state
	 *
	 * @param state state
	 *
	 * return corresponding character
	 */
	public final char getChar(int state)
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
	public final String getTriplet(int state)
	{
		return CODON_TRIPLETS[state];
	}

	/**
	 * Get an array of three nucleotide states making this codon state
	 *
	 * @param state state
	 *
	 * return corresponding triplet string
	 */
	public final int[] getTripletStates(int state)
	{
		int[] triplet = new int[3];
		
		triplet[0] = Nucleotides.INSTANCE.getState(CODON_TRIPLETS[state].charAt(0));
		triplet[1] = Nucleotides.INSTANCE.getState(CODON_TRIPLETS[state].charAt(1));
		triplet[2] = Nucleotides.INSTANCE.getState(CODON_TRIPLETS[state].charAt(2));
		
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
	 * type of data type
	 *
	 * @return integer code for the data type
	 */
	public int getType() {
		return TYPE;
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
     * @return true if the given state represents a stop codon
     */
    public final boolean isStopCodon(int state) {
        return geneticCode.isStopCodon(state);
    }

	// Private members
	
	private GeneticCode geneticCode;

}
