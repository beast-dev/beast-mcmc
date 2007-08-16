/*
 * ConvertAlignment.java
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

package dr.evolution.alignment;

import dr.evolution.datatype.AminoAcids;
import dr.evolution.datatype.CodonTable;
import dr.evolution.datatype.Codons;
import dr.evolution.datatype.DataType;
import dr.evolution.sequence.Sequence;
import dr.evolution.util.Taxon;

/**
 * An alignment class that takes another alignment and converts it on the fly
 * to a different dataType.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 *
 * @version $Id: ConvertAlignment.java,v 1.29 2005/05/24 20:25:55 rambaut Exp $
 */
public class ConvertAlignment extends Alignment.Abstract implements dr.util.XHTMLable
{

	/**
	 * Constructor.
	 */
	public ConvertAlignment() { }

	/**
	 * Constructor.
	 */
	public ConvertAlignment(DataType dataType) {
		this(dataType, null, null);
	}

	/**
	 * Constructor.
	 */
	public ConvertAlignment(DataType dataType, CodonTable codonTable) {
		this(dataType, codonTable, null);
	}

	/**
	 * Constructor.
	 */
	public ConvertAlignment(DataType dataType, Alignment alignment) {
		this(dataType, null, alignment);
	}

	/**
	 * Constructor.
	 */
	public ConvertAlignment(DataType dataType, CodonTable codonTable, Alignment alignment) {
		setDataType(dataType);
		setCodonTable(codonTable);
		setAlignment(alignment);
	}

	/**
	 * Sets the CodonTable of this alignment. 
	 */
	public void setCodonTable(CodonTable codonTable) {
		this.codonTable = codonTable;
	}

	/**
	 * Sets the contained.
	 */
	public void setAlignment(Alignment alignment) {
		if (dataType == null)
			dataType = alignment.getDataType();
			
		this.alignment = alignment;
		
		int newType = dataType.getType();
		int originalType = alignment.getDataType().getType();
		
		if (originalType == DataType.NUCLEOTIDES) {
			if (newType != DataType.CODONS && newType != DataType.AMINO_ACIDS) {
				throw new RuntimeException("Incompatible alignment DataType for ConversionAlignment");
			}
		} else if (originalType == DataType.CODONS) {
			if (newType != DataType.AMINO_ACIDS) {
				throw new RuntimeException("Incompatible alignment DataType for ConversionAlignment");
			}
		} else {
			throw new RuntimeException("Incompatible alignment DataType for ConversionAlignment");
		}		
	}

    // **************************************************************
    // Alignment IMPLEMENTATION
    // **************************************************************

	/**
	 * Sets the dataType of this alignment. This can be different from
	 * the dataTypes of the contained alignment - they will be translated
	 * as required.
	 */
	public void setDataType(DataType dataType) {
		this.dataType = dataType;
	}

	/**
	 * Returns string representation of single sequence in
	 * alignment with gap characters included.
	 */
	public String getAlignedSequenceString(int sequenceIndex) {
		StringBuffer buffer = new StringBuffer();
		for (int i = 0, n = getSiteCount(); i < n; i++) {
			buffer.append(dataType.getChar(getState(sequenceIndex, i)));
		}
		return buffer.toString();
	}


	/**
	 * Returns string representation of single sequence in
	 * alignment with gap characters excluded.
	 */
	public String getUnalignedSequenceString(int sequenceIndex) {
		StringBuffer unaligned = new StringBuffer();
		for (int i = 0, n = getSiteCount(); i < n; i++) {
			
			int state = getState(sequenceIndex, i);
			if (!dataType.isGapState(state)) {
				unaligned.append(dataType.getChar(state));
			}
		}
		
		return unaligned.toString();
	}
	
    // **************************************************************
    // PatternList IMPLEMENTATION
    // **************************************************************

	/**
	 * @return the DataType of this siteList
	 */
	public DataType getDataType() {
		return dataType;
	}
	
    // **************************************************************
    // SiteList IMPLEMENTATION
    // **************************************************************

	/**
	 * @return number of sites
	 */
	public int getSiteCount() {
		if (alignment == null) throw new RuntimeException("ConvertionAlignment has no alignment");
		
		int originalType = alignment.getDataType().getType();
		int count = alignment.getSiteCount();
		
		if (originalType == DataType.NUCLEOTIDES) {
			count /= 3;
		}		
		
		return count;
	}

	/** 
	 * Gets the pattern of site as an array of state numbers (one per sequence) 
	 * @return the site pattern at siteIndex
	 */
	public int[] getSitePattern(int siteIndex) {
		if (alignment == null) throw new RuntimeException("ConvertionAlignment has no alignment");
		
		int i, n = getSequenceCount();
		
		int[] pattern = new int[n];
		
		for (i = 0; i < n; i++) {
			pattern[i] = getState(i, siteIndex);	
		}
		
		return pattern;
	}

	/** 
	 * Gets the pattern index at a particular site
	 * @return the patternIndex
	 */
	public int getPatternIndex(int siteIndex) {
		return siteIndex;
	}

	/** 
	 * @return the sequence state at (taxon, site) 
	 */
	public int getState(int taxonIndex, int siteIndex) {
		if (alignment == null) throw new RuntimeException("ConvertionAlignment has no alignment");

		int newType = dataType.getType();
		int originalType = alignment.getDataType().getType();

		int state = 0;
		
		if (originalType == DataType.NUCLEOTIDES) {
			int siteIndex3 = siteIndex * 3;
			int state1 = alignment.getState(taxonIndex, siteIndex3);
			int state2 = alignment.getState(taxonIndex, siteIndex3 + 1);
			int state3 = alignment.getState(taxonIndex, siteIndex3 + 2);
			
			if (newType == DataType.CODONS) {
				state = Codons.UNIVERSAL.getState(state1, state2, state3);
			} else { // newType == DataType.AMINO_ACIDS
				state = codonTable.getAminoAcidState(Codons.UNIVERSAL.getCanonicalState(Codons.UNIVERSAL.getState(state1, state2, state3)));
			}
			
		} else if (originalType == DataType.CODONS) {
			state = codonTable.getAminoAcidState(alignment.getState(taxonIndex, siteIndex));
		}	
			
		return state;
	}

   // **************************************************************
    // SequenceList IMPLEMENTATION
    // **************************************************************

	/**
	 * @return a count of the number of sequences in the list.
	 */
	public int getSequenceCount() {
		if (alignment == null) throw new RuntimeException("SitePatterns has no alignment");
		return alignment.getSequenceCount();
	}

	/**
	 * @return the ith sequence in the list.
	 */
	public Sequence getSequence(int index) {
		if (alignment == null) throw new RuntimeException("SitePatterns has no alignment");
		return alignment.getSequence(index);
	}

	/**
	 * Sets an named attribute for a given sequence.
	 * @param index the index of the sequence whose attribute is being set.
	 * @param name the name of the attribute.
	 * @param value the new value of the attribute.
	 */
	public void setSequenceAttribute(int index, String name, Object value) {
		if (alignment == null) throw new RuntimeException("SitePatterns has no alignment");
		alignment.setSequenceAttribute(index, name, value);
	}

	/**
	 * @return an object representing the named attributed for the given sequence.
	 * @param index the index of the sequence whose attribute is being fetched.
	 * @param name the name of the attribute of interest.
	 */
	public Object getSequenceAttribute(int index, String name) {
		if (alignment == null) throw new RuntimeException("SitePatterns has no alignment");
		return alignment.getSequenceAttribute(index, name);
	}
	
    // **************************************************************
    // TaxonList IMPLEMENTATION
    // **************************************************************

	/**
	 * @return a count of the number of taxa in the list.
	 */
	public int getTaxonCount() {
		if (alignment == null) throw new RuntimeException("SitePatterns has no alignment");
		return alignment.getTaxonCount();
	}

	/**
	 * @return the ith taxon.
	 */
	public Taxon getTaxon(int taxonIndex) {
		if (alignment == null) throw new RuntimeException("SitePatterns has no alignment");
		return alignment.getTaxon(taxonIndex);
	}

	/**
	 * @return the ID of the ith taxon.
	 */
	public String getTaxonId(int taxonIndex) {
		if (alignment == null) throw new RuntimeException("SitePatterns has no alignment");
		return alignment.getTaxonId(taxonIndex);
	}

	/**
	 * returns the index of the taxon with the given id.
	 */
	public int getTaxonIndex(String id) {
		if (alignment == null) throw new RuntimeException("SitePatterns has no alignment");
		return alignment.getTaxonIndex(id);
	}
	
	/**
	 * returns the index of the given taxon.
	 */
	public int getTaxonIndex(Taxon taxon) {
		if (alignment == null) throw new RuntimeException("SitePatterns has no alignment");
		return alignment.getTaxonIndex(taxon);
	}
	
	/**
	 * @return an object representing the named attributed for the given taxon.
	 * @param taxonIndex the index of the taxon whose attribute is being fetched.
	 * @param name the name of the attribute of interest.
	 */
	public Object getTaxonAttribute(int taxonIndex, String name) {
		if (alignment == null) throw new RuntimeException("SitePatterns has no alignment");
		return alignment.getTaxonAttribute(taxonIndex, name);
	}

	public String toXHTML() {
		String xhtml = "<p><em>Converted Alignment</em> data type = ";
		xhtml += getDataType().getDescription();
		xhtml += ", no. taxa = ";
		xhtml += getTaxonCount();
		xhtml += ", no. sites = ";
		xhtml += getSiteCount();
		xhtml += "</p>";
		
		xhtml += "<pre>";
		
		int length, maxLength = 0;
		for (int i =0; i < getTaxonCount(); i++) {
			length = getTaxonId(i).length();
			if (length > maxLength)
				maxLength = length;
		}
		
		int count, state;
		int type = dataType.getType();
		
		for (int i = 0; i < getTaxonCount(); i++) {
			length = getTaxonId(i).length();
			xhtml += getTaxonId(i);
			for (int j = length; j <= maxLength; j++)
				xhtml += " ";
				
			count = getSiteCount();
			for (int j = 0; j < count; j++) {
				state = getState(i, j);
				if (type == DataType.CODONS)
					xhtml += Codons.UNIVERSAL.getTriplet(state) + " ";
				else
					xhtml += AminoAcids.INSTANCE.getTriplet(state) + " ";
			}
			xhtml += "\n";
		}
		xhtml += "</pre>";
		return xhtml;
	}

	// **************************************************************
	// INSTANCE VARIABLES
	// **************************************************************

	private DataType dataType = null;
	private CodonTable codonTable = null;
	private Alignment alignment = null;
}
