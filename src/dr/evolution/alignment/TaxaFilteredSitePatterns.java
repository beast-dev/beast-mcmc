/*
 * TaxaFilteredSitePatterns.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.evolution.alignment;

import dr.evolution.datatype.DataType;
import dr.evolution.util.Taxon;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

public class TaxaFilteredSitePatterns implements PatternList {

    private PatternList original;
    private String id = null;
    private List<Taxon> taxonList;

    private int[] mapOriginalToFiltered;
    private int[] mapFilteredToOriginal;

    public TaxaFilteredSitePatterns(SitePatterns original, List<Taxon> include, List<Taxon> exclude) {
        this.original = original;

        if (include == null) {
            taxonList = original.asList();
        } else {
            taxonList = include;
        }

        if (exclude != null) {
            for (Taxon t : exclude) {
                if (taxonList.contains(t)) {
                    taxonList.remove(t);
                } else {
                    throw new IllegalArgumentException("Taxon '" + t.getId() + "' is not in current taxon list");
                }
            }
        }

        mapOriginalToFiltered = new int[original.getTaxonCount()];
        mapFilteredToOriginal = new int[this.getTaxonCount()];

        for (int i = 0; i < mapFilteredToOriginal.length; ++i) {
            Taxon taxon = taxonList.get(i);
            mapFilteredToOriginal[i] = original.getTaxonIndex(taxon);
        }

        for (int i = 0; i < mapOriginalToFiltered.length; ++i) {
            Taxon taxon = original.getTaxon(i);
            if (taxonList.contains(taxon)) {
                mapOriginalToFiltered[i] = taxonList.indexOf(taxon);
            } else {
                mapOriginalToFiltered[i] = -1;
            }
        }
    }

    @Override
    public int getPatternCount() {
        return original.getPatternCount();
    }

    @Override
    public int getStateCount() {
        return original.getStateCount();
    }

    @Override
    public int getPatternLength() {
        return original.getPatternLength();
    }

    @Override
    public int[] getPattern(int patternIndex) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double[][] getUncertainPattern(int patternIndex) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public int getPatternState(int taxonIndex, int patternIndex) {

        int originalTaxonIndex = mapFilteredToOriginal[taxonIndex];
        return original.getPatternState(originalTaxonIndex, patternIndex);
    }

    @Override
    public double[] getUncertainPatternState(int taxonIndex, int patternIndex) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double getPatternWeight(int patternIndex) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double[] getPatternWeights() {
        return original.getPatternWeights();

    }

    public int getPatternIndex(int siteIndex){
        // Not implemented yet
        return -1;
    }
    @Override
    public DataType getDataType() {
        return original.getDataType();
    }

    @Override
    public double[] getStateFrequencies() {
        return original.getStateFrequencies();
    }

    @Override
    public boolean areUnique() {
        return original.areUnique();
    }

    @Override
    public boolean areUncertain() {
        return original.areUncertain();
    }

    @Override
    public int getTaxonCount() {
        return taxonList.size();
    }

    @Override
    public Taxon getTaxon(int taxonIndex) {
        return taxonList.get(taxonIndex);
    }

    @Override
    public String getTaxonId(int taxonIndex) {
        return taxonList.get(taxonIndex).getId();
    }

    @Override
    public int getTaxonIndex(String id) {
        for (int i = 0; i < taxonList.size(); ++i) {
            Taxon t = taxonList.get(i);
            if (t.getId().equals(id)) {
                return i;
            }
        }
        return -1;
    }

    @Override
    public int getTaxonIndex(Taxon taxon) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public List<Taxon> asList() {
        List<Taxon> taxa = new ArrayList<Taxon>();
        for (int i = 0, n = getTaxonCount(); i < n; i++) {
            taxa.add(getTaxon(i));
        }
        return taxa;
    }

    @Override
    public Object getTaxonAttribute(int taxonIndex, String name) {
        throw new RuntimeException("Not yet implemented");
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
        throw new RuntimeException("Not yet implemented");
    }
}
