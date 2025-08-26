/*
 * PolymorphismAwareFrequencyModel.java
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
import dr.evolution.datatype.PolymorphismAwareDataType;
import dr.inference.model.Parameter;

/**
 * PoMo frequency model
 *
 * @author Xiang Ji
 * @author Nicola De Maio
 * @author Ben Redelings
 * @author Marc A. Suchard
 */
public class PolymorphismAwareFrequencyModel extends FrequencyModel {

    final PolymorphismAwareDataType dataType;
    final SubstitutionModel baseSubstitution;
    double[] baseQ;

    public PolymorphismAwareFrequencyModel(PolymorphismAwareDataType dataType,
                                           Parameter frequencyParameter,
                                           SubstitutionModel substitutionModel) {
        super(dataType, frequencyParameter);
        this.dataType = dataType;
        this.baseSubstitution = substitutionModel;
        this.baseQ = new double[frequencyParameter.getDimension() * frequencyParameter.getDimension()];
        addModel(substitutionModel);
    }

    public int getFrequencyCount() {
        return dataType.getStateCount();
    }

    private double getWattersonConstant() {
        final int virtualPopSize = dataType.getVirtualPopSize();
        double w = 0;
        for (int i = 1; i < virtualPopSize; i++) {
            w += 1.0 / (double) i;
        }
        return w;
    }

    private double getExpectedNumberMutations() {
        final int baseStateCount = baseSubstitution.getFrequencyModel().getFrequencyCount();

        baseSubstitution.getInfinitesimalMatrix(baseQ);
        double m = 0;
        for (int i = 0; i < baseStateCount; i++) {
            for (int j = 0; j < baseStateCount; j++) {
                if (i != j)
                    m += frequencyParameter.getParameterValue(i) * baseQ[i * baseStateCount + j];
            }
        }
        return m;
    }

    public double[] getFrequencies() {

        final double[] frequencies = new double[dataType.getStateCount()];
        final int virtualPopSize = dataType.getVirtualPopSize();
        final double normalizingConstant = 1.0 / (1.0 + getWattersonConstant() * getExpectedNumberMutations());
        baseSubstitution.getInfinitesimalMatrix(baseQ);
        final int baseStateCount = baseSubstitution.getFrequencyModel().getFrequencyCount();

        for (int i = 0; i < baseStateCount; i++) {
            frequencies[i] = normalizingConstant * frequencyParameter.getParameterValue(i);
        }

        for (int i = 0; i < baseStateCount; i++) {
            for (int j = 0; j < baseStateCount; j++) {
                if (i != j) {
                    for (int k = 1; k < virtualPopSize; k++) {
                        final int stateIndex = dataType.getState(new int[]{i, j}, new int[]{k, virtualPopSize - k});
                        final double q = ((double) k * (virtualPopSize - k)) / (double) virtualPopSize;
                        frequencies[stateIndex] = normalizingConstant * frequencyParameter.getParameterValue(i)
                                * frequencyParameter.getParameterValue(j) * baseQ[i * baseStateCount + j] / q;
                    }
                }
            }
        }

        double sum = 0;
        for (int i = 0; i < frequencies.length; i++) {
            sum += frequencies[i];
        }

        for (int i = 0; i < frequencies.length; i++) {
            frequencies[i] /= sum;
        }

        return frequencies;
    }

}
