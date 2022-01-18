/*
 * TwoParalogGeneConversionSubstitutionModelParser.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
import dr.evolution.datatype.PairedDataType;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;

import java.util.*;

/**
 * @author Xiang Ji
 * @author Jeff Thorne
 * @author Marc A. Suchard
 */
public class PairedParalogSitePatterns implements SiteList {

    private final SitePatterns sitePatterns;
    private final PairedDataType dataType;
    private final Taxa species;
    private final Taxa allSeq;
    private final List<String> paralogs;
    private final String idSeparator;
    private final List<String> singleCopySpecies;
    private final List<String> allIds;

    private int[][] pairedPatterns;
    private String id;

    public PairedParalogSitePatterns(SitePatterns sitePatterns,
                                     String[] paralogs,
                                     String idSeparator,
                                     Taxa species,
                                     String[] singleCopySpecies,
                                     Taxa allSeq) {
        this.sitePatterns = sitePatterns;
        this.dataType = new PairedDataType(sitePatterns.getDataType());
        this.paralogs = Arrays.asList(paralogs);
        this.idSeparator = idSeparator;
        this.species = species;
        this.singleCopySpecies = Arrays.asList(singleCopySpecies);
        this.allSeq = allSeq;
        this.allIds = getAllIds(allSeq);
        this.pairedPatterns = new int[sitePatterns.getPatternCount()][];
        setPairedPatterns();
    }

    private void setPairedPatterns() {

        int patternCount = 0;

        for (int i = 0; i < getSiteCount(); i++) {
            if (sitePatterns.getPatternIndex(i) > patternCount - 1) {
                pairedPatterns[patternCount] = getSitePattern(i);
                patternCount++;
            }
        }
    }

    private List<String> getAllIds(Taxa allSeq) {

        List<String> allIds = new ArrayList<>();

        for (int i = 0; i < allSeq.getTaxonCount(); i++) {
            allIds.add(allSeq.getTaxonId(i));
        }

        return allIds;

    }

    private String getSeqID(String speciesName, int paralogIndex) {
        return speciesName + idSeparator + paralogs.get(paralogIndex);
    }

    @Override
    public int getSiteCount() {
        return sitePatterns.getSiteCount();
    }

    @Override
    public int[] getSitePattern(int siteIndex) {

        int[] pattern = new int[species.getTaxonCount()];

        for (int i = 0; i < species.getTaxonCount(); i++) {
            pattern[i] = getState(i, siteIndex);
        }
        return pattern;
    }

    @Override
    public double[][] getUncertainSitePattern(int siteIndex) {
        throw new UnsupportedOperationException("uncertain patterns not implemented yet");
    }

    @Override
    public int getPatternIndex(int siteIndex) {
        return sitePatterns.getPatternIndex(siteIndex);
    }

    @Override
    public int getState(int taxonIndex, int siteIndex) {
        String speciesName = species.getTaxon(taxonIndex).getId();
        int paralogPattern1, paralogPattern2;

        if (singleCopySpecies.contains(speciesName)) {
            paralogPattern1 = paralogPattern2 = -1;
            for (int j = 0; j < paralogs.size(); j++) {
                final String seqId = getSeqID(speciesName, j);
                if (allIds.contains(seqId)) {
                    final int patternIndex = allSeq.getTaxonIndex(seqId);
                    paralogPattern1 = paralogPattern2 = sitePatterns.getSitePattern(siteIndex)[patternIndex];
                }
            }
            if (paralogPattern1 == -1 && paralogPattern2 == -1) {
                throw new RuntimeException("No matching sequence.");
            }
        } else {
            final String seqId1 = getSeqID(speciesName, 0);
            final String seqId2 = getSeqID(speciesName, 1);
            final int patternIndex1 = allSeq.getTaxonIndex(seqId1);
            final int patternIndex2 = allSeq.getTaxonIndex(seqId2);
            paralogPattern1 = sitePatterns.getSitePattern(siteIndex)[patternIndex1];
            paralogPattern2 = sitePatterns.getSitePattern(siteIndex)[patternIndex2];
        }

        return dataType.getState(paralogPattern1, paralogPattern2);
    }

    @Override
    public double[] getUncertainState(int taxonIndex, int siteIndex) {
        throw new UnsupportedOperationException("uncertain patterns not implemented yet");
    }

    @Override
    public int getPatternCount() {
        return sitePatterns.getPatternCount();
    }

    @Override
    public int getStateCount() {
        return dataType.getStateCount();
    }

    @Override
    public int getPatternLength() {
        return species.getTaxonCount();
    }

    @Override
    public int[] getPattern(int patternIndex) {
        return pairedPatterns[patternIndex];
    }

    @Override
    public double[][] getUncertainPattern(int patternIndex) {
        throw new UnsupportedOperationException("uncertain patterns not implemented yet");
    }

    @Override
    public int getPatternState(int taxonIndex, int patternIndex) {
        return pairedPatterns[patternIndex][taxonIndex];
    }

    @Override
    public double[] getUncertainPatternState(int taxonIndex, int patternIndex) {
        throw new UnsupportedOperationException("uncertain patterns not implemented yet");
    }

    @Override
    public double getPatternWeight(int patternIndex) {
        return sitePatterns.getPatternWeight(patternIndex);
    }

    @Override
    public double[] getPatternWeights() {
        return sitePatterns.getPatternWeights();
    }

    @Override
    public DataType getDataType() {
        return dataType;
    }

    @Override
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

    @Override
    public int getTaxonCount() {
        return species.getTaxonCount();
    }

    @Override
    public Taxon getTaxon(int taxonIndex) {
        return species.getTaxon(taxonIndex);
    }

    @Override
    public String getTaxonId(int taxonIndex) {
        return species.getTaxonId(taxonIndex);
    }

    @Override
    public int getTaxonIndex(String id) {
        return species.getTaxonIndex(id);
    }

    @Override
    public int getTaxonIndex(Taxon taxon) {
        return species.getTaxonIndex(taxon);
    }

    @Override
    public List<Taxon> asList() {
        return species.asList();
    }

    @Override
    public Object getTaxonAttribute(int taxonIndex, String name) {
        return species.getTaxonAttribute(taxonIndex, name);
    }

    @Override
    public String getId() {
        return id;
    }

    @Override
    public void setId(String id) {
        this.id = id;
    }

    @Override
    public Iterator<Taxon> iterator() {
        return species.iterator();
    }
}
