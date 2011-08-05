/*
 * SitePatterns.java
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
import dr.evolution.datatype.Codons;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;

import java.util.*;

/**
 * Stores a set of site patterns. This differs from the simple Patterns
 * class because it stores the pattern index for each site. Thus it has
 * a connection to a single alignment.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: SitePatterns.java,v 1.47 2005/06/22 16:44:17 rambaut Exp $
 */
public class SitePatterns implements SiteList, dr.util.XHTMLable {

    /**
     * the source alignment
     */
    protected SiteList siteList = null;

    /**
     * number of sites
     */
    protected int siteCount = 0;

    /**
     * number of patterns
     */
    protected int patternCount = 0;

    /**
     * length of site patterns
     */
    protected int patternLength = 0;

    /**
     * site -> site pattern
     */
    protected int[] sitePatternIndices;

    /**
     * count of invariant patterns
     */
    protected int invariantCount;

    /**
     * weights of each site pattern
     */
    protected double[] weights;

    /**
     * site patterns [site pattern][taxon]
     */
    protected int[][] patterns;

    protected int from, to, every;

    protected boolean strip = true;  // Strip out completely ambiguous sites

    protected boolean unique = true; // Compress into weighted list of unique patterns

    /**
     * Constructor
     */
    public SitePatterns(Alignment alignment) {
        this(alignment, null, 0, 0, 1);
    }

    /**
     * Constructor
     */
    public SitePatterns(Alignment alignment, TaxonList taxa) {
        this(alignment, taxa, 0, 0, 1);
    }

    /**
     * Constructor
     */
    public SitePatterns(Alignment alignment, int from, int to, int every) {
        this(alignment, null, from, to, every);
    }

//    /**
//     * Constructor for dnds
//     */
//    public SitePatterns(Alignment alignment, int from, int to, int every, boolean unique) {
//        this(alignment, null, from, to, every, unique);
//    }
    
    /**
     * Constructor
     */

    public SitePatterns(Alignment alignment, TaxonList taxa, int from, int to, int every) {
        this(alignment,taxa,from,to,every,true);
    }

    public SitePatterns(Alignment alignment, TaxonList taxa, int from, int to, int every, boolean strip) {
        this(alignment, taxa, from, to, every, strip, true);
    }

    public SitePatterns(Alignment alignment, TaxonList taxa, int from, int to, int every, boolean strip, boolean unique) {
        if (taxa != null) {
            SimpleAlignment a = new SimpleAlignment();

            for (int i = 0; i < alignment.getSequenceCount(); i++) {
                if (taxa.getTaxonIndex(alignment.getTaxonId(i)) != -1) {
                    a.addSequence(alignment.getSequence(i));
                }
            }

            alignment = a;
        }
        this.strip = strip;
        this.unique = unique;
        
        setPatterns(alignment, from, to, every);
    }

    /**
     * Constructor
     */
    public SitePatterns(SiteList siteList) {
        this(siteList, -1, -1, 1);
    }

    /**
     * Constructor
     */
    public SitePatterns(SiteList siteList, int from, int to, int every) {
        setPatterns(siteList, from, to, every);
    }

    public SiteList getSiteList() {
        return siteList;
    }

    public int getFrom() {
        return from;
    }

    public int getTo() {
        return to;
    }

    public int getEvery() {
        return every;
    }

    public void setFrom(int from) {
        setPatterns(getSiteList(), from, getTo(), getEvery());
    }

    public void setTo(int to) {
        setPatterns(getSiteList(), getFrom(), to, getEvery());
    }

    public void setEvery(int every) {
        setPatterns(getSiteList(), getFrom(), getTo(), every);
    }

    /**
     * sets up pattern list using an alignment
     */
    public void setPatterns(SiteList siteList, int from, int to, int every) {

        this.siteList = siteList;
        this.from = from;
        this.to = to;
        this.every = every;

        if (siteList == null) {
            return;
        }

        if (from <= -1)
            from = 0;

        if (to <= -1)
            to = siteList.getSiteCount() - 1;

        if (every <= 0)
            every = 1;

        siteCount = ((to - from) / every) + 1;

        patternCount = 0;

        patterns = new int[siteCount][];

        sitePatternIndices = new int[siteCount];
        weights = new double[siteCount];

        invariantCount = 0;
        int[] pattern;

        int site = 0;

        for (int i = from; i <= to; i += every) {
            pattern = siteList.getSitePattern(i);

            if (!strip || !isInvariant(pattern) ||
                    (!isGapped(pattern) &&
                            !isAmbiguous(pattern) &&
                            !isUnknown(pattern))) {

                sitePatternIndices[site] = addPattern(pattern);

            }  else {
              sitePatternIndices[site] = -1;
            }
            site++;
        }
    }

    /**
     * adds a pattern to the pattern list
     *
     * @return the index of the pattern in the pattern list
     */
    private int addPattern(int[] pattern) {

        for (int i = 0; i < patternCount; i++) {

            if (unique && comparePatterns(patterns[i], pattern)) {

                weights[i] += 1.0;
                return i;
            }
        }

        if (isInvariant(pattern)) {
            invariantCount++;
        }

        int index = patternCount;
        patterns[index] = pattern;
        weights[index] = 1.0;
        patternCount++;

        return index;
    }

    /**
     * @return true if the pattern is invariant
     */
    private boolean isGapped(int[] pattern) {
        int len = pattern.length;

        for (int i = 0; i < len; i++) {
            if (getDataType().isGapState(pattern[i])) {
                return true;
            }
        }
        return false;
    }

    /**
     * @return true if the pattern is invariant
     */
    private boolean isAmbiguous(int[] pattern) {
        int len = pattern.length;

        for (int i = 0; i < len; i++) {
            if (getDataType().isAmbiguousState(pattern[i])) {
                return true;
            }
        }
        return false;
    }

    /**
     * @return true if the pattern is invariant
     */
    private boolean isUnknown(int[] pattern) {
        int len = pattern.length;

        for (int i = 0; i < len; i++) {
            if (getDataType().isUnknownState(pattern[i])) {
                return true;
            }
        }
        return false;
    }

    /**
     * @return true if the pattern is invariant
     */
    private boolean isInvariant(int[] pattern) {
        int len = pattern.length;

        int state = pattern[0];
        for (int i = 1; i < len; i++) {
            if (pattern[i] != state) {
                return false;
            }
        }

        return true;
    }

    /**
     * compares two patterns
     *
     * @return true if they are identical
     */
    protected boolean comparePatterns(int[] pattern1, int[] pattern2) {

        int len = pattern1.length;
        for (int i = 0; i < len; i++) {
            if (pattern1[i] != pattern2[i]) {
                return false;
            }
        }

        return true;
    }

    /**
     * @return number of invariant sites (these will be first in the list).
     */
    public int getInvariantCount() {
        return invariantCount;
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
     *
     * @return the site pattern at siteIndex
     */
    public int[] getSitePattern(int siteIndex) {
        final int sitePatternIndice = sitePatternIndices[siteIndex];
        return sitePatternIndice >= 0 ? patterns[sitePatternIndice] : null;
    }

    /**
     * Gets the pattern index at a particular site
     *
     * @return the patternIndex
     */
    public int getPatternIndex(int siteIndex) {
        return sitePatternIndices[siteIndex];
    }

    /**
     * @return the sequence state at (taxon, site)
     */
    public int getState(int taxonIndex, int siteIndex) {
        final int sitePatternIndice = sitePatternIndices[siteIndex];
        // is that right?
        return sitePatternIndice >= 0 ? patterns[sitePatternIndice][taxonIndex] : getDataType().getGapState();
    }

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
        if (siteList == null) throw new RuntimeException("SitePatterns has no alignment");
        return siteList.getStateCount();
    }

    /**
     * Gets the length of the pattern strings which will usually be the
     * same as the number of taxa
     *
     * @return the length of patterns
     */
    public int getPatternLength() {
        return getTaxonCount();
    }

    /**
     * Gets the pattern as an array of state numbers (one per sequence)
     *
     * @return the pattern at patternIndex
     */
    public int[] getPattern(int patternIndex) {
        return patterns[patternIndex];
    }

    /**
     * @return state at (taxonIndex, patternIndex)
     */
    public int getPatternState(int taxonIndex, int patternIndex) {
        return patterns[patternIndex][taxonIndex];
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
        if (siteList == null) throw new RuntimeException("SitePatterns has no alignment");
        return siteList.getDataType();
    }

    /**
     * @return the frequency of each state
     */
    public double[] getStateFrequencies() {
        return PatternList.Utils.empiricalStateFrequencies(this);
    }

    // **************************************************************
    // TaxonList IMPLEMENTATION
    // **************************************************************

    /**
     * @return a count of the number of taxa in the list.
     */
    public int getTaxonCount() {
        if (siteList == null) throw new RuntimeException("SitePatterns has no alignment");
        return siteList.getTaxonCount();
    }

    /**
     * @return the ith taxon.
     */
    public Taxon getTaxon(int taxonIndex) {
        if (siteList == null) throw new RuntimeException("SitePatterns has no alignment");
        return siteList.getTaxon(taxonIndex);
    }

    /**
     * @return the ID of the ith taxon.
     */
    public String getTaxonId(int taxonIndex) {
        if (siteList == null) throw new RuntimeException("SitePatterns has no alignment");
        return siteList.getTaxonId(taxonIndex);
    }

    /**
     * returns the index of the taxon with the given id.
     */
    public int getTaxonIndex(String id) {
        if (siteList == null) throw new RuntimeException("SitePatterns has no alignment");
        return siteList.getTaxonIndex(id);
    }

    /**
     * returns the index of the given taxon.
     */
    public int getTaxonIndex(Taxon taxon) {
        if (siteList == null) throw new RuntimeException("SitePatterns has no alignment");
        return siteList.getTaxonIndex(taxon);
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

    /**
     * @param taxonIndex the index of the taxon whose attribute is being fetched.
     * @param name       the name of the attribute of interest.
     * @return an object representing the named attributed for the given taxon.
     */
    public Object getTaxonAttribute(int taxonIndex, String name) {
        if (siteList == null) throw new RuntimeException("SitePatterns has no alignment");
        return siteList.getTaxonAttribute(taxonIndex, name);
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
        String xhtml = "<p><em>Pattern List</em>  pattern count = ";
        xhtml += getPatternCount();
        xhtml += "  invariant count = ";
        xhtml += getInvariantCount();
        xhtml += "</p>";

        xhtml += "<pre>";

        int count, state;
        int type = getDataType().getType();

        count = getPatternCount();

        int length, maxLength = 0;
        for (int i = 0; i < count; i++) {
            length = Integer.toString((int) getPatternWeight(i)).length();
            if (length > maxLength)
                maxLength = length;
        }

        for (int i = 0; i < count; i++) {
            length = Integer.toString(i + 1).length();
            for (int j = length; j < maxLength; j++)
                xhtml += " ";
            xhtml += Integer.toString(i + 1) + ": ";

            length = Integer.toString((int) getPatternWeight(i)).length();
            xhtml += Integer.toString((int) getPatternWeight(i));
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
