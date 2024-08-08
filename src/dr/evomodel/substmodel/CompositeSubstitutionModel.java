/*
 * GeneralSubstitutionModel.java
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

package dr.evomodel.substmodel;

import dr.evolution.datatype.DataType;
import dr.inference.model.DuplicatedParameter;
import dr.inference.model.Parameter;

import java.util.ArrayList;
import java.util.List;

/**
 * A substitution model that is a composite of two other models.
 *
 * @author Andrew Rambaut
 */
public class CompositeSubstitutionModel extends BaseSubstitutionModel {

    /**
     * the rate which the others are set relative to
     */
    protected int ratesRelativeTo;

    public CompositeSubstitutionModel(String name, DataType dataType,
                                      SubstitutionModel substitutionModel1,
                                      SubstitutionModel substitutionModel2) {
        super(name, dataType, substitutionModel1.getFrequencyModel(), null);

        this.substitutionModel1 = substitutionModel1;
        this.substitutionModel2 = substitutionModel2;

        addModel(substitutionModel1);
        addModel(substitutionModel2);
    }

    protected void frequenciesChanged() {
        // Nothing to precalculate
    }

    protected void ratesChanged() {
        // Nothing to precalculate
    }

    protected void setupRelativeRates(double[] rates) {
        for (int i = 0; i < rates.length; i++) {
            if (i == ratesRelativeTo) {
                rates[i] = 1.0;
            } else if (ratesRelativeTo < 0 || i < ratesRelativeTo) {
                rates[i] = ratesParameter.getParameterValue(i);
            } else {
                rates[i] = ratesParameter.getParameterValue(i - 1);
            }
        }
    }

    public void setNormalization(boolean doNormalization) {
        this.doNormalization = doNormalization;
    }

    protected double getNormalizationValue(double[][] matrix, double[] pi) {
        double norm = 1.0;
        if (doNormalization) {
            norm = super.getNormalizationValue(matrix, pi);
        }
        return norm;
    }


    // *****************************************************************
    // Interface Model
    // *****************************************************************

    protected void storeState() {
    } // nothing to do

    /**
     * Restore the additional stored state
     */
    protected void restoreState() {
        updateMatrix = true;
    }

    protected void acceptState() {
    } // nothing to do

    private final SubstitutionModel substitutionModel1;
    private final SubstitutionModel substitutionModel2;


    protected Parameter ratesParameter = null;
    private boolean doNormalization = true;
}