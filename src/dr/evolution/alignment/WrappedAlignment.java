/*
 * WrappedAlignment.java
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
import dr.evolution.sequence.Sequence;
import dr.evolution.util.Taxon;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id$
 */
public abstract class WrappedAlignment extends Alignment.Abstract {
    protected WrappedAlignment(Alignment alignment) {
        this.alignment = alignment;
    }

    /**
     * @return the sequence state at (taxon, site)
     */
    public abstract int getState(int taxonIndex, int siteIndex);

    @Override
    public double[] getUncertainState(int taxonIndex, int siteIndex) {
        throw new UnsupportedOperationException("getUncertainState not implemented yet");
    }


    public void setDataType(DataType dataType) {
        // do nothing by default
    }

    /**
     * Returns string representation of single sequence in
     * alignment with gap characters included.
     */
    public String getAlignedSequenceString(int sequenceIndex) {
        StringBuffer buffer = new StringBuffer();
        for (int i = 0, n = getSiteCount(); i < n; i++) {
            buffer.append(getDataType().getChar(getState(sequenceIndex, i)));
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
            if (!getDataType().isGapState(state)) {
                unaligned.append(getDataType().getChar(state));
            }
        }

        return unaligned.toString();
    }

    @Override
    public double[][] getUncertainPattern(int patternIndex) {
        throw new UnsupportedOperationException("getUncertainPattern not implemented yet");
    }

    @Override
    public double[] getUncertainPatternState(int taxonIndex, int patternIndex) {
        throw new UnsupportedOperationException("getUncertainPatternState not implemented yet");
    }

    /**
     * @return the DataType of this siteList
     */
    public DataType getDataType() {
        return alignment.getDataType();
    }

    @Override
    public boolean areUncertain() {
        return false;
    }

    /**
     * @return number of sites
     */
    public int getSiteCount() {
        return alignment.getSiteCount();
    }

    /**
     * Gets the pattern of site as an array of state numbers (one per sequence)
     * @return the site pattern at siteIndex
     */
    public int[] getSitePattern(int siteIndex) {
        int i, n = getSequenceCount();

        int[] pattern = new int[n];

        for (i = 0; i < n; i++) {
            pattern[i] = getState(i, siteIndex);
        }

        return pattern;
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
        return alignment.getPatternIndex(siteIndex);
    }

    /**
     * @return a count of the number of sequences in the list.
     */
    public int getSequenceCount() {
        return alignment.getSequenceCount();
    }

    /**
     * @return the ith sequence in the list.
     */
    public Sequence getSequence(int index) {
        return alignment.getSequence(index);
    }

    /**
     * Sets an named attribute for a given sequence.
     * @param index the index of the sequence whose attribute is being set.
     * @param name the name of the attribute.
     * @param value the new value of the attribute.
     */
    public void setSequenceAttribute(int index, String name, Object value) {
        alignment.setSequenceAttribute(index, name, value);
    }

    /**
     * @return an object representing the named attributed for the given sequence.
     * @param index the index of the sequence whose attribute is being fetched.
     * @param name the name of the attribute of interest.
     */
    public Object getSequenceAttribute(int index, String name) {
        return alignment.getSequenceAttribute(index, name);
    }

    /**
     * @return a count of the number of taxa in the list.
     */
    public int getTaxonCount() {
        return alignment.getTaxonCount();
    }

    /**
     * @return the ith taxon.
     */
    public Taxon getTaxon(int taxonIndex) {
        return alignment.getTaxon(taxonIndex);
    }

    /**
     * @return the ID of the ith taxon.
     */
    public String getTaxonId(int taxonIndex) {
        return alignment.getTaxonId(taxonIndex);
    }

    /**
     * returns the index of the taxon with the given id.
     */
    public int getTaxonIndex(String id) {
        return alignment.getTaxonIndex(id);
    }

    /**
     * returns the index of the given taxon.
     */
    public int getTaxonIndex(Taxon taxon) {
        return alignment.getTaxonIndex(taxon);
    }

    public List<Taxon> asList() {
        List<Taxon> taxa = new ArrayList<Taxon>();
        for (int i = 0, n = getTaxonCount(); i < n; i++) {
            taxa.add(getTaxon(i));
        }
        return taxa;
    }

    public String toString() {
        dr.util.NumberFormatter formatter = new dr.util.NumberFormatter(6);

        StringBuffer buffer = new StringBuffer();

        for (int i = 0; i < getSequenceCount(); i++) {
            String name = formatter.formatToFieldWidth(getTaxonId(i), 10);
            buffer.append(">").append(name).append("\n");
            buffer.append(getAlignedSequenceString(i)).append("\n");
        }

        return buffer.toString();
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
     * @return an object representing the named attributed for the given taxon.
     * @param taxonIndex the index of the taxon whose attribute is being fetched.
     * @param name the name of the attribute of interest.
     */
    public Object getTaxonAttribute(int taxonIndex, String name) {
        return alignment.getTaxonAttribute(taxonIndex, name);
    }

    protected Alignment alignment = null;

}
