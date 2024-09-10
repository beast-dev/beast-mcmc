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
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

/**
 * A substitution model that is a composite of two other models.
 *
 * @author Andrew Rambaut
 */
public class CompositeSubstitutionModel extends ComplexSubstitutionModel {

    /**
     *
     * @param name
     * @param dataType
     * @param substitutionModel1
     * @param weightParameter1 the weight (rate) of the first substitution model
     * @param substitutionModel2
     * @param weightParameter2 the weight (rate) of the second substitution model
     */
    public CompositeSubstitutionModel(String name, DataType dataType,
                                      GeneralSubstitutionModel substitutionModel1,
                                      Parameter weightParameter1,
                                      GeneralSubstitutionModel substitutionModel2,
                                      Parameter weightParameter2
    ) {
        super(name, dataType, substitutionModel1.getFrequencyModel(), null);

        this.substitutionModel1 = substitutionModel1;
        int rateCount1 = substitutionModel1.rateCount;
        this.substitutionModel2 = substitutionModel2;
        this.weightParameter1 = weightParameter1;
        int rateCount2 = substitutionModel2.rateCount;
        this.weightParameter2 = weightParameter2;

        this.rateCount = Math.max(rateCount1, rateCount2);
        addModel(substitutionModel1);
        addModel(substitutionModel2);
        addVariable(weightParameter1);
        addVariable(weightParameter2);

        setStateCount(dataType.getStateCount());
    }

    @Override
    protected int getRateCount(int stateCount) {
        if (substitutionModel1 == null) {
            return super.getRateCount(stateCount);
        }
        int rateCount1 = substitutionModel1.rateCount;
        int rateCount2 = substitutionModel2.rateCount;
        return Math.max(rateCount1, rateCount2);
    }

    protected void frequenciesChanged() {
        // Nothing to precalculate
    }

    protected void ratesChanged() {
        // Nothing to precalculate
    }

    @Override
    protected void setupRelativeRates(double[] rates) {
        double weight1 = weightParameter1.getParameterValue(0);
        double weight2 = weightParameter2.getParameterValue(0);
        double[] rates1 = substitutionModel1.getRelativeRates();
        double[] rates2 = substitutionModel2.getRelativeRates();
        if (substitutionModel1.rateCount == substitutionModel2.rateCount) {
            // both 6 or 12 rate
            for (int i = 0; i < rateCount; i++) {
                rates[i] = (rates1[i] * weight1) + (rates2[i] * weight2);
            }
        } else if (substitutionModel1.rateCount < substitutionModel2.rateCount) {
            // substitutionModel1 is 6 rate, substitutionModel2 12 rate
            int k = 0;
            for (int i = 0; i < substitutionModel1.rateCount; i++) {
                rates[k] = (rates1[i] * weight1) + (rates2[k] * weight2);
                k++;
            }
            for (int i = 0; i < substitutionModel1.rateCount; i++) {
                rates[k] = (rates1[i] * weight1) + (rates2[k] * weight2);
                k++;
            }
        } else {
            // substitutionModel1 is 12 rate, substitutionModel2 is 6 rate
            int k = 0;
            for (int i = 0; i < substitutionModel2.rateCount; i++) {
                rates[k] = (rates1[k] * weight1) + (rates2[i] * weight2);
                k++;
            }
            for (int i = 0; i < substitutionModel2.rateCount; i++) {
                rates[k] = (rates1[k] * weight1) + (rates2[i] * weight2);
                k++;
            }
        }
    }

    protected void setupQMatrix(double[] rates, double[] pi, double[][] matrix) {

        int i, j, k = 0;
        for (i = 0; i < stateCount; i++) {
            for (j = i + 1; j < stateCount; j++) {
                matrix[i][j] = rates[k] * pi[j];
                k++;
            }
        }

        if (k == rates.length) {
            // this is a symmetrical matrix (6 rates) so restart the rate indexing
            k = 0;
        }
        // Copy lower triangle in column-order form (transposed)
        for (j = 0; j < stateCount; j++) {
            for (i = j + 1; i < stateCount; i++) {
                matrix[i][j] = rates[k] * pi[j];
                k++;
            }
        }
    }

    // *****************************************************************
    // Interface Model
    // *****************************************************************

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        updateMatrix = true;
        fireModelChanged();
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        updateMatrix = true;
        fireModelChanged();
    }

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

    private final GeneralSubstitutionModel substitutionModel1;
    private final GeneralSubstitutionModel substitutionModel2;


    protected final Parameter weightParameter1;
    protected final Parameter weightParameter2;
    private boolean doNormalization = true;
}