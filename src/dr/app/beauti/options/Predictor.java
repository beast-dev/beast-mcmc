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
import java.util.*;

/**
 * @author Andrew Rambaut
 */
public class Predictor implements Serializable {
    private static final long serialVersionUID = -9152518508699327745L;

    public enum PredictorType {
        MATRIX,
        FROM_VECTOR,
        TO_VECTOR,
        BOTH_VECTOR;

        public String toString() {
            return name().toLowerCase();
        }
    }

    private final PredictorType predictorType;
    private final TraitData trait;

    private String name;
    private boolean isIncluded;
    private boolean isLogged;
    private boolean isStandardized;

    protected final BeautiOptions options;

    public Predictor(BeautiOptions options, String name, TraitData trait, Map<String, Double> vector, PredictorType predictorType) {
        this.options = options;
        this.name = name;
        this.trait = trait;
        this.predictorType = predictorType;
    }

    public Predictor(BeautiOptions options, String name, TraitData trait, Map<String, List<Double>> matrix) {
        this.options = options;
        this.name = name;
        this.trait = trait;
        this.predictorType = PredictorType.MATRIX;
    }

    /////////////////////////////////////////////////////////////////////////

    public PredictorType getType() {
        return predictorType;
    }

    public String getName() {
        return name;
    }

    public void setName(final String name) {
        this.name = name;
    }

    public boolean isIncluded() {
        return isIncluded;
    }

    public void setIncluded(boolean included) {
        isIncluded = included;
    }

    public boolean isLogged() {
        return isLogged;
    }

    public void setLogged(boolean logged) {
        isLogged = logged;
    }

    public boolean isStandardized() {
        return isStandardized;
    }

    public void setStandardized(boolean standardized) {
        isStandardized = standardized;
    }

    public String toString() {
        return name;
    }
}
