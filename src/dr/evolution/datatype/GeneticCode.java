/*
 * GeneticCode.java
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
 * A set of standard genetic codes.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 *
 * @version $Id: GeneticCode.java,v 1.11 2005/05/24 20:25:56 rambaut Exp $
 */

public final class GeneticCode implements CodonTable {

	public static final String GENETIC_CODE = "geneticCode";
	
	/**
	 * Constants used to refer to the built in code tables
	 */
	public static final int UNIVERSAL_ID = 0;
	public static final int VERTEBRATE_MT_ID = 1;
	public static final int YEAST_ID = 2;
	public static final int MOLD_PROTOZOAN_MT_ID = 3;
	public static final int MYCOPLASMA_ID = 4;
	public static final int INVERTEBRATE_MT_ID = 5;
	public static final int CILIATE_ID = 6;
	public static final int ECHINODERM_MT_ID = 7;
	public static final int EUPLOTID_NUC_ID = 8;
	public static final int BACTERIAL_ID = 9;
	public static final int ALT_YEAST_ID = 10;
	public static final int ASCIDIAN_MT_ID = 11;
	public static final int FLATWORM_MT_ID = 12;
	public static final int BLEPHARISMA_NUC_ID = 13;
	public static final int NO_STOPS_ID = 14;

	/**
	 * Standard genetic code tables from GENBANK
	 * Nucleotides go A, C, G, T - Note: this is not the order used by the Genbank web site
	 * With the first codon position most significant (i.e. AAA, AAC, AAG, AAT, ACA, etc.).
	 */
	public static final String[] GENETIC_CODE_TABLES = {
		// Universal
		"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF",
		// Vertebrate Mitochondrial
		"KNKNTTTT*S*SMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF",
		// Yeast
		"KNKNTTTTRSRSMIMIQHQHPPPPRRRRTTTTEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF",
		// Mold Protozoan Mitochondrial
		"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF",
		// Mycoplasma
		"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF",
		// Invertebrate Mitochondrial
		"KNKNTTTTSSSSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF",
		// Ciliate
		"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVQYQYSSSS*CWCLFLF",
		// Echinoderm Mitochondrial
		"NNKNTTTTSSSSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF",
		// Euplotid Nuclear
		"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSCCWCLFLF",
		// Bacterial
		"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF",
		// Alternative Yeast
		"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLSLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF",
		// Ascidian Mitochondrial
		"KNKNTTTTGSGSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF",
		// Flatworm Mitochondrial
		"NNKNTTTTSSSSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVYY*YSSSSWCWCLFLF",
		// Blepharisma Nuclear
		"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*YQYSSSS*CWCLFLF",
		// No stops
		"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVYYQYSSSSWCWCLFLF"
	};

	/**
	 * Names of the standard genetic code tables from GENBANK
	 */
	public static final String[] GENETIC_CODE_NAMES = {
		"universal", "vertebrateMitochondrial", "yeast", "moldProtozoanMitochondrial",
		"mycoplasma", "invertebrateMitochondrial", "ciliate", "echinodermMitochondrial",
		"euplotidNuclear", "bacterial", "alternativeYeast", "ascidianMitochondrial",
		"flatwormMitochondrial", "blepharismaNuclear", "noStops"
	};

	/**
	 * Descriptions of the standard genetic code tables from GENBANK
	 */
	public static final String[] GENETIC_CODE_DESCRIPTIONS = {
		"Universal", "Vertebrate Mitochondrial", "Yeast", "Mold Protozoan Mitochondrial",
		"Mycoplasma", "Invertebrate Mitochondrial", "Ciliate", "Echinoderm Mitochondrial",
		"Euplotid Nuclear", "Bacterial", "Alternative Yeast", "Ascidian Mitochondrial",
		"Flatworm Mitochondrial", "Blepharisma Nuclear", "Test case with no stop codons"
	};

	public static final GeneticCode UNIVERSAL = new GeneticCode(UNIVERSAL_ID);
	public static final GeneticCode VERTEBRATE_MT = new GeneticCode(VERTEBRATE_MT_ID);
	public static final GeneticCode YEAST = new GeneticCode(YEAST_ID);
	public static final GeneticCode MOLD_PROTOZOAN_MT = new GeneticCode(MOLD_PROTOZOAN_MT_ID);
	public static final GeneticCode MYCOPLASMA = new GeneticCode(MYCOPLASMA_ID);
	public static final GeneticCode INVERTEBRATE_MT = new GeneticCode(INVERTEBRATE_MT_ID);
	public static final GeneticCode CILIATE = new GeneticCode(CILIATE_ID);
	public static final GeneticCode ECHINODERM_MT = new GeneticCode(ECHINODERM_MT_ID);
	public static final GeneticCode EUPLOTID_NUC = new GeneticCode(EUPLOTID_NUC_ID);
	public static final GeneticCode BACTERIAL = new GeneticCode(BACTERIAL_ID);
	public static final GeneticCode ALT_YEAST = new GeneticCode(ALT_YEAST_ID);
	public static final GeneticCode ASCIDIAN_MT = new GeneticCode(ASCIDIAN_MT_ID);
	public static final GeneticCode FLATWORM_MT = new GeneticCode(FLATWORM_MT_ID);
	public static final GeneticCode BLEPHARISMA_NUC = new GeneticCode(BLEPHARISMA_NUC_ID);
	public static final GeneticCode NO_STOPS = new GeneticCode(NO_STOPS_ID);

	public static final GeneticCode[] GENETIC_CODES = {
		UNIVERSAL, VERTEBRATE_MT, YEAST, MOLD_PROTOZOAN_MT, MYCOPLASMA, INVERTEBRATE_MT,
		CILIATE, ECHINODERM_MT, EUPLOTID_NUC, BACTERIAL, ALT_YEAST, ASCIDIAN_MT,
		FLATWORM_MT, BLEPHARISMA_NUC, NO_STOPS
	};

	public GeneticCode(int geneticCode) {
		
		this.geneticCode = geneticCode;
		codeTable = GENETIC_CODE_TABLES[geneticCode];
	}
	
	/**
	 * Returns the name of the genetic code
	 */
	public String getName() {
		return GENETIC_CODE_NAMES[geneticCode];
	}
	
	/**
	 * Returns the description of the genetic code
	 */
	public String getDescription() {
		return GENETIC_CODE_DESCRIPTIONS[geneticCode];
	}
	
	/**
	 * Returns the char associated with AminoAcid represented by codonState.
	 * Note that the char is as defined by AminoAcids.java
	 * @see AminoAcids
	 * @see Codons
	 * @return state for '?' if codon unknown
	 */
	public char getAminoAcidChar(int codonState) {
		if (codonState == Codons.UNKNOWN_STATE)
			return AminoAcids.UNKNOWN_CHARACTER;
		else if (codonState == Codons.GAP_STATE)
			return AminoAcids.GAP_CHARACTER;
			
		return codeTable.charAt(codonState);
	}
	
	/**
	 * Returns the state associated with AminoAcid represented by codonState.
	 * Note that the state is the canonical state (generated combinatorially)
	 * @see AminoAcids
	 * @see Codons
	 * @return '?' if codon unknown
	 */
	public int getAminoAcidState(int codonState) {
		if (codonState == Codons.UNKNOWN_STATE)
			return AminoAcids.UNKNOWN_STATE;
		else if (codonState == Codons.GAP_STATE)
			return AminoAcids.GAP_STATE;
			
		return AminoAcids.AMINOACID_STATES[getAminoAcidChar(codonState)];
	}

	/**
	 * Note that the state is the canonical state (generated combinatorially)
     * @return whether the codonState is a stop codon
	 */
	public boolean isStopCodon(int codonState) {
		return (getAminoAcidState(codonState) == AminoAcids.STOP_STATE);
	}

	/**
	 * @return all the possible codons for a given amino acid
	 */
	public char[][] getCodonsFromAminoAcidState(int aminoAcidState) {
		throw new RuntimeException("not yet implemented");
	}

	/*
	 * @return all the possible codons for a given amino acid
	 */
	public char[][] getCodonsFromAminoAcidChar(char aminoAcidChar) {
		throw new RuntimeException("not yet implemented");
	}

	/*
	 * @returns three IUPAC states representing the given amino acid
	 * @note The returned array should not be altered, and implementations
	 *       should attempt to implement this as efficiently as possible
	 */
	public int[] getAmbiguousCodonFromAminoAcidState(int aminoAcid) {
		throw new RuntimeException("not yet implemented");
	}

	/**
	 * @return the codon states of stop amino acids.
	 */
	public int[] getStopCodonIndices() {
	
		int i, j, n = getStopCodonCount();
		int[] indices = new int[n];
		
		j = 0;
		for (i = 0; i < 64; i++) {
			if (codeTable.charAt(i) == AminoAcids.STOP_CHARACTER) {
				indices[j] = i;
				j++;
			}
		}
		
		return indices;
	}

	/**
	 * Returns the number of terminator amino acids.
	 */
	public int getStopCodonCount() {
		int i, count = 0;
		
		for (i = 0; i < 64; i++) {
			if (codeTable.charAt(i) == AminoAcids.STOP_CHARACTER)
				count++;
		}
		
		return count;
	}
	
	private int geneticCode;
	private String codeTable;

}
