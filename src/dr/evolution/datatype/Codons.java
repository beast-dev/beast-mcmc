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
 * @version $Id: Codons.java,v 1.15 2005/05/24 20:25:56 rambaut Exp $
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 */
public class Codons extends DataType {
    
	/**
	 * Name of data type. For XML and human reading of data type.
	 */
	public static final String DESCRIPTION = "codon";
	public static final int TYPE = CODONS;
	
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
		
		stateCount = 64 - geneticCode.getStopCodonCount();
		ambiguousStateCount = 66;
		
		stateMap = new int[ambiguousStateCount];
		reverseMap = new int[ambiguousStateCount];
		
		int j = 0;
		int k = stateCount;
		for (int i = 0; i < 64; i++) {
			if (!geneticCode.isStopCodon(i)) {
				stateMap[j] = i;
				reverseMap[i] = j;
				j++;
			} else {
				stateMap[k] = i;
				reverseMap[i] = k;
				k++;
			}
		}
		for (int i = 64; i < ambiguousStateCount; i++) {
			stateMap[i] = i;
			reverseMap[i] = i;
		}
		
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
	 * @return the canonical state (in standard combinatorial order)
	 * of a funny codon state.
	 */
	public final int getCanonicalState(int funnyState) {
		return stateMap[funnyState];

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
	
		int canonicalState = (nuc1 * 16) + (nuc2 * 4) + nuc3;

		return reverseMap[canonicalState];
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
		return CODON_TRIPLETS[stateMap[state]];
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
		
		triplet[0] = Nucleotides.INSTANCE.getState(CODON_TRIPLETS[stateMap[state]].charAt(0));
		triplet[1] = Nucleotides.INSTANCE.getState(CODON_TRIPLETS[stateMap[state]].charAt(1));
		triplet[2] = Nucleotides.INSTANCE.getState(CODON_TRIPLETS[stateMap[state]].charAt(2));
		
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
     * @return true if the given state is a stop codon
     */
    public final boolean isStopCodon(int state) {
        return geneticCode.isStopCodon(stateMap[state]);
    }

	// Private members
	
	private GeneticCode geneticCode;
	private int[] stateMap;
	private int[] reverseMap;

    static {
        registerDataType(Codons.DESCRIPTION+"-"+GeneticCode.UNIVERSAL.getName(), Codons.UNIVERSAL);
        registerDataType(Codons.DESCRIPTION+"-"+GeneticCode.VERTEBRATE_MT.getName(), Codons.VERTEBRATE_MT);
        registerDataType(Codons.DESCRIPTION+"-"+GeneticCode.YEAST.getName(), Codons.YEAST);
        registerDataType(Codons.DESCRIPTION+"-"+GeneticCode.MOLD_PROTOZOAN_MT.getName(), Codons.MOLD_PROTOZOAN_MT);
        registerDataType(Codons.DESCRIPTION+"-"+GeneticCode.MYCOPLASMA.getName(), Codons.MYCOPLASMA);
        registerDataType(Codons.DESCRIPTION+"-"+GeneticCode.INVERTEBRATE_MT.getName(), Codons.INVERTEBRATE_MT);
        registerDataType(Codons.DESCRIPTION+"-"+GeneticCode.CILIATE.getName(), Codons.CILIATE);
        registerDataType(Codons.DESCRIPTION+"-"+GeneticCode.ECHINODERM_MT.getName(), Codons.ECHINODERM_MT);
        registerDataType(Codons.DESCRIPTION+"-"+GeneticCode.EUPLOTID_NUC.getName(), Codons.EUPLOTID_NUC);
        registerDataType(Codons.DESCRIPTION+"-"+GeneticCode.BACTERIAL.getName(), Codons.BACTERIAL);
        registerDataType(Codons.DESCRIPTION+"-"+GeneticCode.ALT_YEAST.getName(), Codons.ALT_YEAST);
        registerDataType(Codons.DESCRIPTION+"-"+GeneticCode.ASCIDIAN_MT.getName(), Codons.ASCIDIAN_MT);
        registerDataType(Codons.DESCRIPTION+"-"+GeneticCode.FLATWORM_MT.getName(), Codons.FLATWORM_MT);
        registerDataType(Codons.DESCRIPTION+"-"+GeneticCode.BLEPHARISMA_NUC.getName(), Codons.BLEPHARISMA_NUC);
    }
}
