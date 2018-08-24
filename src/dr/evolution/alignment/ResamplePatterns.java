/*
 * ResamplePatterns.java
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

package dr.evolution.alignment;

import dr.evolution.datatype.AminoAcids;
import dr.evolution.datatype.Codons;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.util.Taxon;

import java.util.*;

/**
 * Provides an abstract base class for resampling patterns (i.e. bootstrap and jackknife).
 *
 * @version $Id: ResamplePatterns.java,v 1.12 2005/05/24 20:25:55 rambaut Exp $
 *
 * @author Andrew Rambaut
 */
public abstract class ResamplePatterns implements PatternList, dr.util.XHTMLable
{
	/** the source patterns */
	protected SiteList patterns = null;

	/** number of site patterns */
	protected int patternCount = 0;

	/** weights of each site pattern */
	protected double[] weights;

	/** bootstrapped pattern -> source pattern */
	protected int[] patternIndices;

	/**
	 * Sets the source of patterns
	 */
	public void setPatterns(SiteList patterns) {
		this.patterns = patterns;
	}

	/**
	 * Perform a resampling of the patterns
	 */
	public abstract void resamplePatterns();

    // **************************************************************
    // PatternList IMPLEMENTATION
    // **************************************************************

	/**
	 * @return number of patterns
	 */
	public int getPatternCount() {
		return patternCount;
	}

	/**
	 * @return number of states for this siteList
	 */
	public int getStateCount() {
		if (patterns == null) throw new RuntimeException("ResamplePatterns has no source patterns");
		return patterns.getStateCount();
	}

	/**
	 * Gets the length of the pattern strings which will usually be the
	 * same as the number of taxa
	 * @return the length of patterns
	 */
	public int getPatternLength() {
		if (patterns == null) throw new RuntimeException("ResamplePatterns has no source patterns");
		return patterns.getPatternLength();
	}

	/**
	 * Gets the pattern as an array of state numbers (one per sequence)
	 * @return the pattern at patternIndex
	 */
	public int[] getPattern(int patternIndex) {
		if (patterns == null) throw new RuntimeException("ResamplePatterns has no source patterns");
		return patterns.getPattern(patternIndices[patternIndex]);
	}

	@Override
	public double[][] getUncertainPattern(int patternIndex) {
		if (patterns == null) throw new RuntimeException("ResamplePatterns has no source patterns");
		return patterns.getUncertainPattern(patternIndices[patternIndex]);
	}

	/**
	 * @return state at (taxonIndex, patternIndex)
	 */
	public int getPatternState(int taxonIndex, int patternIndex) {
		return getPattern(patternIndex)[taxonIndex];
	}

	@Override
	public double[] getUncertainPatternState(int taxonIndex, int patternIndex) {
		return getUncertainPattern(patternIndex)[taxonIndex];
	}

	/**
	 * Gets the weight of a site pattern
	 */
	public double getPatternWeight(int patternIndex) {
		return weights[patternIndex];
	}

	/**
	 * @return the array of pattern weights
	 */
	public double[] getPatternWeights() {
		return weights;
	}

	/**
	 * @return the DataType of this siteList
	 */
	public DataType getDataType() {
		if (patterns == null) throw new RuntimeException("ResamplePatterns has no source patterns");
		return patterns.getDataType();
	}

	/**
	 * @return the frequency of each state
	 */
	public double[] getStateFrequencies() {
		return PatternList.Utils.empiricalStateFrequencies(this);
	}

	@Override
	public boolean areUnique() {
		return patterns.areUnique();
	}

	@Override
	public boolean areUncertain() {
		return patterns.areUncertain();
	}

	// **************************************************************
    // TaxonList IMPLEMENTATION
    // **************************************************************

	/**
	 * @return a count of the number of taxa in the list.
	 */
	public int getTaxonCount() {
		if (patterns == null) throw new RuntimeException("ResamplePatterns has no source patterns");
		return patterns.getTaxonCount();
	}

	/**
	 * @return the ith taxon.
	 */
	public Taxon getTaxon(int taxonIndex) {
		if (patterns == null) throw new RuntimeException("ResamplePatterns has no source patterns");
		return patterns.getTaxon(taxonIndex);
	}

	/**
	 * @return the ID of the ith taxon.
	 */
	public String getTaxonId(int taxonIndex) {
		if (patterns == null) throw new RuntimeException("ResamplePatterns has no source patterns");
		return patterns.getTaxonId(taxonIndex);
	}

	/**
	 * returns the index of the taxon with the given id.
	 */
	public int getTaxonIndex(String id) {
		if (patterns == null) throw new RuntimeException("ResamplePatterns has no source patterns");
		return patterns.getTaxonIndex(id);
	}

	/**
	 * returns the index of the given taxon.
	 */
	public int getTaxonIndex(Taxon taxon) {
		if (patterns == null) throw new RuntimeException("SitePatterns has no alignment");
		return patterns.getTaxonIndex(taxon);
	}

	/**
	 * @return an object representing the named attributed for the given taxon.
	 * @param taxonIndex the index of the taxon whose attribute is being fetched.
	 * @param name the name of the attribute of interest.
	 */
	public Object getTaxonAttribute(int taxonIndex, String name) {
		if (patterns == null) throw new RuntimeException("ResamplePatterns has no source patterns");
		return patterns.getTaxonAttribute(taxonIndex, name);
	}

    public List<Taxon> asList() {
        List<Taxon> taxa = new ArrayList<Taxon>();
        for (int i = 0, n = getTaxonCount(); i < n; i++) {
            taxa.add(getTaxon(i));
        }
        return taxa;
    }

    public Iterator<Taxon> iterator() {
        return new Iterator<Taxon>() {
            private int index = -1;

            public boolean hasNext() {
                return index < getTaxonCount() - 1;
            }

            public Taxon next() {
                index ++;
                return getTaxon(index);
            }

            public void remove() { /* do nothing */ }
        };
    }
    // **************************************************************
    // Identifiable IMPLEMENTATION
    // **************************************************************

	protected String id = null;

	/**
	 * @return the id.
	 */
	public String getId() {
		return id;
	}

	/**
	 * Sets the id.
	 */
	public void setId(String id) {
		this.id = id;
	}

    // **************************************************************
    // XHTMLable IMPLEMENTATION
    // **************************************************************

	public String toXHTML() {
		String xhtml = "<p><em>Jackknife Pattern List</em>  pattern count = ";
		xhtml += getPatternCount();
		xhtml += "</p>";

		xhtml += "<pre>";

		int count, state;
		int type = getDataType().getType();

		count = getPatternCount();

		int length, maxLength = 0;
		for (int i = 0; i < count; i++) {
			length = Integer.toString((int)getPatternWeight(i)).length();
			if (length > maxLength)
				maxLength = length;
		}

		for (int i = 0; i < count; i++) {
			length = Integer.toString(i+1).length();
			for (int j = length; j < maxLength; j++)
				xhtml += " ";
			xhtml += Integer.toString(i+1) + ": ";

			length = Integer.toString((int)getPatternWeight(i)).length();
			xhtml += Integer.toString((int)getPatternWeight(i));
			for (int j = length; j <= maxLength; j++)
				xhtml += " ";

			for (int j = 0; j < getTaxonCount(); j++) {
				state = getPatternState(j, i);

				if (type == DataType.NUCLEOTIDES) {
					xhtml += Nucleotides.INSTANCE.getChar(state) + " ";
				} else if (type == DataType.CODONS) {
					xhtml += Codons.UNIVERSAL.getTriplet(state) + " ";
				} else {
					xhtml += AminoAcids.INSTANCE.getChar(state) + " ";
				}
			}
			xhtml += "\n";
		}
		xhtml += "</pre>";
		return xhtml;
	}

}
