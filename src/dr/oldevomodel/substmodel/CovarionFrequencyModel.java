/*
 * CovarionFrequencyModel.java
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

package dr.oldevomodel.substmodel;

import dr.evolution.datatype.DataType;
import dr.inference.model.Parameter;

/**
 * @author Alexei Drummond
 */
public class CovarionFrequencyModel extends FrequencyModel {

    public CovarionFrequencyModel(DataType dataType, Parameter frequencyParameter, Parameter hiddenFrequencies) {
        super(dataType, frequencyParameter);

        this.hiddenFrequencies = hiddenFrequencies;
        addVariable(hiddenFrequencies);
    }

    public double[] getFrequencies() {

        int k = 0;

        int numStates = frequencyParameter.getDimension();
        int numHiddenStates = hiddenFrequencies.getDimension();

        double[] freqs = new double[numStates * numHiddenStates];
        for (int i = 0; i < numHiddenStates; i++) {
            for (int j = 0; j < numStates; j++) {
                freqs[k] = frequencyParameter.getParameterValue(j) *
                        hiddenFrequencies.getParameterValue(i);
                k += 1;
            }
        }

        return freqs;
    }

    public void setFrequency(int i, double value) {
        throw new UnsupportedOperationException();
    }

    public double getFrequency(int i) {

        int numStates = frequencyParameter.getDimension();

        return frequencyParameter.getParameterValue(i % numStates) *
                hiddenFrequencies.getParameterValue(i / numStates);
    }

    public int getFrequencyCount() {
        return frequencyParameter.getDimension() * hiddenFrequencies.getDimension();
    }

    Parameter hiddenFrequencies;
}
