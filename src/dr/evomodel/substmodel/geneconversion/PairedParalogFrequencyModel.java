/*
 * PairedParalogFrequencyModel.java
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

package dr.evomodel.substmodel.geneconversion;

import dr.evolution.datatype.PairedDataType;
import dr.evomodel.substmodel.FrequencyModel;
import dr.inference.model.Parameter;

/**
 * @author Xiang Ji
 * @author Jeff Thorne
 * @author Marc A. Suchard
 */
public class PairedParalogFrequencyModel extends FrequencyModel {

    PairedDataType dataType;

    public PairedParalogFrequencyModel(PairedDataType dataType, Parameter frequencyParameter) {

        super(dataType, frequencyParameter);

        this.dataType = dataType;
    }

    public double[] getFrequencies() {

        final double[] frequencies = new double[dataType.getStateCount()];

        for (int i = 0; i < getFrequencyParameter().getDimension(); i++) {
            frequencies[dataType.getState(i, i)] = getFrequencyParameter().getParameterValue(i);
        }

        return frequencies;
    }

    public void setFrequency(int i, double value) {
        // do nothing
    }

    public int getFrequencyCount() {
        return dataType.getStateCount();
    }

    public double getFrequency(int i) {
        final int numStates = getFrequencyParameter().getDimension();

        final int index1 = getState1(i, numStates);
        final int index2 = getState2(i, numStates);

        double frequency = 0.0;

        if (index1 == index2) {
            frequency = getFrequencyParameter().getParameterValue(index1);
        }

        return frequency;
    }

    public int getState1(int state, int numState) {
        return state / numState;
    }


    public int getState2(int state, int numState) {
        return state % numState;
    }
}
