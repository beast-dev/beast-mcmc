/*
 * TraitData.java
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

package dr.app.beauti.options;

import dr.evolution.datatype.ContinuousDataType;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.GeneralDataType;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @author Walter Xie
 */
public class TraitData implements Serializable {
    private static final long serialVersionUID = -9152518508699327745L;

    public enum TraitType {
        DISCRETE,
        INTEGER,
        CONTINUOUS;

        public String toString() {
            return name().toLowerCase();
        }
    }

    private TraitType traitType = TraitType.DISCRETE;

    private final String fileName;
    private String name;

    protected final BeautiOptions options;

    private List<Predictor> predictorList = new ArrayList<Predictor>();

    public TraitData(BeautiOptions options, String name, String fileName, TraitType traitType) {
        this.options = options;
        this.name = name;
        this.fileName = fileName;
        this.traitType = traitType;
    }

    /////////////////////////////////////////////////////////////////////////

    public TraitType getTraitType() {
        return traitType;
    }

    public void setTraitType(TraitType traitType) {
        this.traitType = traitType;
    }

//    public TraitOptions getTraitOptions() {
//        return traitOptions;
//    }

    public int getSiteCount() {
        return 0;
    }

    public int getTaxaCount() {
        return options.taxonList.getTaxonCount();
    }

    public Taxon getTaxon(int i) {
        return options.taxonList.getTaxon(i);
    }

    public boolean hasValue(int i) {
        if (options.taxonList == null || options.taxonList.getTaxon(i) == null
                || options.taxonList.getTaxon(i).getAttribute(getName()) == null) return false;
        return options.taxonList.getTaxon(i).getAttribute(getName()).toString().trim().length() > 0;
    }

    public DataType
    getDataType() {
        switch (traitType) {
            case DISCRETE:
                return GeneralDataType.INSTANCE;
            case INTEGER:
                return GeneralDataType.INSTANCE;
            case CONTINUOUS:
                return ContinuousDataType.INSTANCE;
        }
        throw new IllegalArgumentException("Unknown trait type");
    }

    public String getFileName() {
        return fileName;
    }

    public String getName() {
        return name;
    }

    public void setName(final String name) {
        this.name = name;
    }

    public Set<String> getStatesOfTrait() {
        return getStatesListOfTrait(options.taxonList, getName());
    }

    public Set<String> getStatesOfTrait(Taxa taxonList) {
        return getStatesListOfTrait(taxonList, getName());
    }


    public static Set<String> getStatesListOfTrait(Taxa taxonList, String traitName) {
        Set<String> states = new TreeSet<String>();

        if (taxonList == null) {
            throw new IllegalArgumentException("taxon list is null");
        }

        for (int i = 0; i < taxonList.getTaxonCount(); i++) {
            Taxon taxon = taxonList.getTaxon(i);
            String attr = (String) taxon.getAttribute(traitName);

            // ? is used to denote missing data so is not a state...
            if (attr != null && !attr.equals("?")) {
                states.add(attr);
            }
        }


        return states;
    }

    public static int getEmptyStateIndex(Taxa taxonList, String traitName) {
        if (taxonList == null) {
            throw new IllegalArgumentException("taxon list is null");
        }

        for (int i = 0; i < taxonList.getTaxonCount(); i++) {
            Taxon taxon = taxonList.getTaxon(i);
            String attr = (String) taxon.getAttribute(traitName);

            // ? is used to denote missing data so is not a state...
            if (attr == null || attr.equals("?")) {
                return i;
            }
        }

        return -1;
    }

    public void addPredictor(Predictor predictor) {
        predictorList.add(predictor);
    }

    public void removePredictor(Predictor predictor) {
        predictorList.remove(predictor);
    }

    /**
     * Returns the number of included predictors (counting origin and destination ones separately)
     * @return
     */
    public int getIncludedPredictorCount() {
        int count = 0;
        for (Predictor predictor : getPredictors()) {
            if (predictor.isIncluded()) {
                count +=  1 + (predictor.getType() == Predictor.Type.BOTH_VECTOR ? 1 : 0);
            }
        }
        return count;
    }

    public List<Predictor> getIncludedPredictors() {
        List<Predictor> includedPredictors = new ArrayList<Predictor>();
        for (Predictor predictor : getPredictors()) {
            if (predictor.isIncluded()) {
                includedPredictors.add(predictor);
            }
        }
        return includedPredictors;
    }

    public List<Predictor> getPredictors() {
        return predictorList;
    }

    public String toString() {
        return name;
    }
}
