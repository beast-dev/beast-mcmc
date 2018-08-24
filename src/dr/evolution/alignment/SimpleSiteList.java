/*
 * SimpleSiteList.java
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

import dr.evolution.datatype.*;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;

import java.util.List;
import java.util.Iterator;

/**
 * Stores a set of site patterns. This is setup by loading site patterns using the
 * addPattern method. It is assumed they are added in order.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 *
 * @version $Id: SimpleSiteList.java,v 1.1 2005/06/23 13:53:41 rambaut Exp $
 */
public class SimpleSiteList implements SiteList {

	private final TaxonList taxonList;
	private final DataType dataType;

    private int siteCount = 0;
    private int[][] sitePatterns = new int[0][];

	public SimpleSiteList(DataType dataType) {
        this.taxonList = null;
        this.dataType = dataType;
    }

	public SimpleSiteList(DataType dataType, TaxonList taxonList) {
		this.taxonList = taxonList;
        this.dataType = dataType;
	}

    public int addPattern(int[] pattern) {
        int capacity = sitePatterns.length;
        if (siteCount >= capacity) {
            capacity += 10000;
            int[][] newSitePatterns = new int[capacity][];
            for (int i = 0; i < siteCount; i++) {
                newSitePatterns[i] = sitePatterns[i];
            }
            sitePatterns = newSitePatterns;
        }

        sitePatterns[siteCount] = pattern;
        siteCount++;

        return siteCount - 1;
    }

	// **************************************************************
	// SiteList IMPLEMENTATION
	// **************************************************************

	/**
	 * @return number of sites
	 */
	public int getSiteCount() {
		return siteCount;
	}

	/**
	 * Gets the pattern of site as an array of state numbers (one per sequence)
	 * @return the site pattern at siteIndex
	 */
	public int[] getSitePattern(int siteIndex) {
		return sitePatterns[siteIndex];
	}

	@Override
	public double[][] getUncertainSitePattern(int siteIndex) {
		throw new UnsupportedOperationException("getUncertainSitePattern not implemented yet");
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
		return getSitePattern(siteIndex)[taxonIndex];
	}

	@Override
	public double[] getUncertainState(int taxonIndex, int siteIndex) {
		throw new UnsupportedOperationException("getUncertainState not implemented yet");
	}

	// **************************************************************
	// PatternList IMPLEMENTATION
	// **************************************************************

	/**
	 * @return number of patterns
	 */
	public int getPatternCount() {
		return getSiteCount();
	}

	/**
	 * @return number of states for this siteList
	 */
	public int getStateCount() {
		return dataType.getStateCount();
	}

	/**
	 * Gets the length of the pattern strings which will usually be the
	 * same as the number of taxa
	 * @return the length of patterns
	 */
	public int getPatternLength() {
		return getTaxonCount();
	}

	/**
	 * Gets the pattern as an array of state numbers (one per sequence)
	 * @return the pattern at patternIndex
	 */
	public int[] getPattern(int patternIndex) {
		return getSitePattern(patternIndex);
	}

	@Override
	public double[][] getUncertainPattern(int patternIndex) {
		throw new UnsupportedOperationException("getUncertainPattern not implemented yet");
	}

	/**
	 * @return state at (taxonIndex, patternIndex)
	 */
	public int getPatternState(int taxonIndex, int patternIndex) {
		return getSitePattern(patternIndex)[taxonIndex];
	}

	@Override
	public double[] getUncertainPatternState(int taxonIndex, int patternIndex) {
		throw new UnsupportedOperationException("getUncertainPatternState not implemented yet");
	}

	/**
	 * Gets the weight of a site pattern
	 */
	public double getPatternWeight(int patternIndex) {
		return 1.0;
	}

	/**
	 * @return the array of pattern weights
	 */
	public double[] getPatternWeights() {
		double[] weights = new double[getSiteCount()];
		for (int i = 0; i < getSiteCount(); i++) {
			weights[i] = 1.0;
		}
		return weights;
	}

	/**
	 * @return the DataType of this siteList
	 */
	public DataType getDataType() {
		return dataType;
	}

	/**
	 * @return the frequency of each state
	 */
	public double[] getStateFrequencies() {
		return Utils.empiricalStateFrequencies(this);
	}

	@Override
	public boolean areUnique() {
		return false;
	}

	@Override
	public boolean areUncertain() {
		return false;
	}

	// **************************************************************
	// TaxonList IMPLEMENTATION
	// **************************************************************

	/**
	 * @return a count of the number of taxa in the list.
	 */
	public int getTaxonCount() {
		if (taxonList == null) throw new RuntimeException("SimpleSiteList has no taxonList");
		return taxonList.getTaxonCount();
	}

	/**
	 * @return the ith taxon.
	 */
	public Taxon getTaxon(int taxonIndex) {
		if (taxonList == null) throw new RuntimeException("SimpleSiteList has no taxonList");
		return taxonList.getTaxon(taxonIndex);
	}

	/**
	 * @return the ID of the ith taxon.
	 */
	public String getTaxonId(int taxonIndex) {
		if (taxonList == null) throw new RuntimeException("SimpleSiteList has no taxonList");
		return taxonList.getTaxonId(taxonIndex);
	}

	/**
	 * returns the index of the taxon with the given id.
	 */
	public int getTaxonIndex(String id) {
		if (taxonList == null) throw new RuntimeException("SimpleSiteList has no taxonList");
		return taxonList.getTaxonIndex(id);
	}

	/**
	 * returns the index of the given taxon.
	 */
	public int getTaxonIndex(Taxon taxon) {
		if (taxonList == null) throw new RuntimeException("SimpleSiteList has no taxonList");
		return taxonList.getTaxonIndex(taxon);
	}

    public List<Taxon> asList() {
        if (taxonList == null) throw new RuntimeException("SimpleSiteList has no taxonList");
        return taxonList.asList();
    }

    public Iterator<Taxon> iterator() {
        if (taxonList == null) throw new RuntimeException("SimpleSiteList has no taxonList");
        return taxonList.iterator();
    }


    /**
	 * @return an object representing the named attributed for the given taxon.
	 * @param taxonIndex the index of the taxon whose attribute is being fetched.
	 * @param name the name of the attribute of interest.
	 */
	public Object getTaxonAttribute(int taxonIndex, String name) {
		if (taxonList == null) throw new RuntimeException("SimpleSiteList has no taxonList");
		return taxonList.getTaxonAttribute(taxonIndex, name);
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
}
