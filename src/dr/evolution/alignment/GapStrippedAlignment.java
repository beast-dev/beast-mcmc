/*
 * GapStrippedAlignment.java
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

import dr.evolution.datatype.DataType;
import dr.evolution.sequence.Sequence;
import dr.evolution.util.Taxon;

import java.util.*;

/**
 * Provides bootstrap replicate patterns
 *
 * @version $Id: GapStrippedAlignment.java,v 1.4 2005/05/24 20:25:55 rambaut Exp $
 *
 * @author Andrew Rambaut
 */
public class GapStrippedAlignment extends Alignment.Abstract
{
    Alignment alignment;
    boolean[] hasGap = null;

    /**
     * Constructor
     */
    public GapStrippedAlignment(Alignment a) {
        this.alignment = a;
        DataType dataType = a.getDataType();
        hasGap = new boolean[a.getSiteCount()];
        for (int i = 0; i < hasGap.length; i++) {
            for (int j = 0; j < a.getSequenceCount(); j++) {
                if (dataType.isGapState(alignment.getState(j,i))) {
                    hasGap[i] = true;
                    break;
                }
            }
        }
    }

    public final void setDataType(DataType dataType) {
        throw new UnsupportedOperationException();
        //alignment.setDataType(dataType);
    }

    public final String getAlignedSequenceString(int sequenceIndex) {
        return getSequence(sequenceIndex).getSequenceString();
    }

    public final String getUnalignedSequenceString(int sequenceIndex) {
        return getSequence(sequenceIndex).getSequenceString();
    }

    public final int getSequenceCount() {
        return alignment.getSequenceCount();
    }

    /**
     * Very inefficient implementation, use sparingly
     * @param sequenceIndex
     * @return
     */
    public final Sequence getSequence(int sequenceIndex) {

        DataType dataType = getDataType();

        StringBuffer buffer = new StringBuffer();
        for (int i = 0; i < alignment.getSiteCount(); i++) {
            if (!hasGap[i]) {
                buffer.append(dataType.getChar(alignment.getState(sequenceIndex,i)));
            }
        }
        return new Sequence(buffer.toString());
    }

    public final void setSequenceAttribute(int index, String name, Object value) {
        throw new UnsupportedOperationException();
    }

    public final Object getSequenceAttribute(int index, String name) {
        throw new UnsupportedOperationException();
    }

    public final int getTaxonCount() {
        return alignment.getTaxonCount();
    }

    public final Taxon getTaxon(int taxonIndex) {
        return alignment.getTaxon(taxonIndex);
    }

    public final String getTaxonId(int taxonIndex) {
        return alignment.getTaxonId(taxonIndex);
    }

    public final int getTaxonIndex(String id) {
        return alignment.getTaxonIndex(id);
    }

    public final int getTaxonIndex(Taxon taxon) {
        return alignment.getTaxonIndex(taxon);
    }

    public final Object getTaxonAttribute(int taxonIndex, String name) {
        return alignment.getTaxonAttribute(taxonIndex, name);
    }

    public List<Taxon> asList() {
        return alignment.asList();
    }

    public Iterator<Taxon> iterator() {
        return alignment.iterator();
    }

    public final int getSiteCount() {
        int siteCount = 0;
        for (int i = 0; i < hasGap.length; i++) {
            if (!hasGap[i]) siteCount += 1;
        }
        return siteCount;
    }

    public final int[] getSitePattern(int siteIndex) {
        return alignment.getSitePattern(fullIndex(siteIndex));
    }

    @Override
    public double[][] getUncertainSitePattern(int siteIndex) {
        return alignment.getUncertainSitePattern(fullIndex(siteIndex));
    }

    public final int getPatternIndex(int siteIndex) {
        return alignment.getPatternIndex(fullIndex(siteIndex));
    }

    public final int getState(int taxonIndex, int siteIndex) {

        return alignment.getState(taxonIndex, fullIndex(siteIndex));
    }

    @Override
    public double[] getUncertainState(int taxonIndex, int siteIndex) {
        return alignment.getUncertainState(taxonIndex, fullIndex(siteIndex));
    }

    @Override
    public double[][] getUncertainPattern(int patternIndex) {
        return alignment.getUncertainPattern(fullIndex(patternIndex));
    }

    @Override
    public double[] getUncertainPatternState(int taxonIndex, int patternIndex) {
        return new double[0];
    }

    public final DataType getDataType() {
        return alignment.getDataType();
    }

    @Override
    public boolean areUncertain() {
        return false;
    }

    private final int fullIndex(int gapStrippedIndex) {
        int index = 0;
        int fullIndex = 0;
        while (index < gapStrippedIndex) {
            if (!hasGap[fullIndex]) index += 1;
            fullIndex += 1;
        }
        return fullIndex;
    }
}
