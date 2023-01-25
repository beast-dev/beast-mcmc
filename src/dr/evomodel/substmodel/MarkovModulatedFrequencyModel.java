/*
 * MarkovModulatedFrequencyModel.java
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

package dr.evomodel.substmodel;

import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.LUDecomposition;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

import java.util.List;

/**
 * @author Marc A. Suchard
 * @author Vladimir Minin
 *         <p/>
 *         A class for implementing a kronecker sum of CTMC models in BEAST using BEAGLE
 *         This work is supported by NSF grant 0856099
 *         <p/>
 *         O'Brien JD, Minin VN and Suchard MA (2009) Learning to count: robust estimates for labeled distances between
 *         molecular sequences. Molecular Biology and Evolution, 26, 801-814
 */

public class MarkovModulatedFrequencyModel extends FrequencyModel {

    MarkovModulatedFrequencyModel(String name, List<FrequencyModel> freqModels, Parameter switchingRates,
                                  Parameter relativeWeight) {
        super(name);
        this.freqModels = freqModels;
        int freqCount = 0;
        stateCount = freqModels.get(0).getFrequencyCount();
        numBaseModel = freqModels.size();
        for (int i = 0; i < numBaseModel; i++) {
            int size = freqModels.get(i).getFrequencyCount();
            if (stateCount != size) {
                throw new RuntimeException("MarkovModulatedFrequencyModel requires all frequencies model to have the same dimension");
            }
            addModel(freqModels.get(i));

            freqCount += size;
        }
        totalFreqCount = freqCount;
        this.switchingRates = switchingRates;
        addVariable(switchingRates);
        
        baseStationaryDistribution = new double[numBaseModel];
        storedBaseStationaryDistribution = new double[numBaseModel];
        stationaryDistributionKnown = false;

        DoubleMatrix2D d = new DenseDoubleMatrix2D(numBaseModel, numBaseModel);
        d.set(0, 0, 1.0);

        this.relativeWeight = relativeWeight;
        if (relativeWeight != null) {
            addVariable(relativeWeight);
        }
        checkRelativeWeight();
    }

    public void setFrequency(int i, double value) {
        throw new RuntimeException("Not implemented");
    }

    public double getFrequency(int index) {
        int whichModel = index / stateCount;
        int whichState = index % stateCount;
        double relativeFreq = freqModels.get(whichModel).getFrequency(whichState);

        if (relativeWeight != null) {
            relativeFreq *= relativeWeight.getParameterValue(whichModel);
        } else {
            relativeFreq /= numBaseModel;
        }

        // Scale by stationary distribution over hidden classes
        if (numBaseModel > 1) {
            if (!stationaryDistributionKnown) {
                computeStationaryDistribution(baseStationaryDistribution);
                stationaryDistributionKnown = true;
            }
        }

        return relativeFreq;
    }

    private void computeStationaryDistribution(double[] stationaryDistribution) {

        if (allRatesAreZero(switchingRates)) {
            return;
        }

        // Uses an LU decomposition to solve Q^t \pi = 0 and \sum \pi_i = 1
        DoubleMatrix2D mat2 = new DenseDoubleMatrix2D(numBaseModel + 1, numBaseModel);
        int index2 = 0;
        for (int i = 0; i < numBaseModel; ++i) {
            for (int j = i + 1; j < numBaseModel; ++j) {
                mat2.set(j, i, switchingRates.getParameterValue(index2)); // Transposed
                index2++;
            }
        }
        for (int j = 0; j < numBaseModel; ++j) {
            for (int i = j + 1; i < numBaseModel; ++i) {
                mat2.set(j, i, switchingRates.getParameterValue(index2)); // Transposed
                index2++;
            }
        }
        for (int i = 0; i < numBaseModel; ++i) {
            double rowTotal = 0.0;
            for (int j = 0; j < numBaseModel; ++j) {
                if (i != j) {
                    rowTotal += mat2.get(j, i); // Transposed
                }

            }
            mat2.set(i, i, -rowTotal);
        }

        // Add row for sum-to-one constraint
        for (int i = 0; i < numBaseModel; ++i) {
            mat2.set(numBaseModel, i, 1.0);
        }

        LUDecomposition decomposition = new LUDecomposition(mat2);
        DoubleMatrix2D x = new DenseDoubleMatrix2D(numBaseModel + 1, 1);
        x.set(numBaseModel, 0, 1.0);
        DoubleMatrix2D y = decomposition.solve(x);
        for (int i = 0; i < numBaseModel; ++i) {
            stationaryDistribution[i] = y.get(i, 0);
        }
    }

    private static boolean allRatesAreZero(Parameter rates) {
        for (int i = 0; i < rates.getDimension(); ++i) {
            if (rates.getParameterValue(i) != 0.0) {
                return false;
            }
        }
        return true;
    }

    protected void storeState() {
        System.arraycopy(baseStationaryDistribution, 0, storedBaseStationaryDistribution, 0, numBaseModel);
        storedStationaryDistributionKnown = stationaryDistributionKnown;
    }

    protected void restoreState() {
        double[] tmp = baseStationaryDistribution;
        baseStationaryDistribution = storedBaseStationaryDistribution;
        storedBaseStationaryDistribution = tmp;

        stationaryDistributionKnown = storedStationaryDistributionKnown;
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if (variable == switchingRates) {
            stationaryDistributionKnown = false;
        }
    }

    private void checkRelativeWeight() {
        if (relativeWeight != null) {
            double sum = 0.0;
            for (double x : relativeWeight.getParameterValues()) {
                sum += x;
            }
            if (sum != 1.0) {
                throw new IllegalArgumentException("Relative weights must sum to 1.0");
            }
        }
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        fireModelChanged();
    }

    public int getFrequencyCount() {
        return totalFreqCount;
    }

    public Parameter getFrequencyParameter() {
        throw new RuntimeException("Not implemented");
    }

    private List<FrequencyModel> freqModels;

    private final int numBaseModel;
    private final int totalFreqCount;
    private final int stateCount;
    private final Parameter switchingRates;
    private final Parameter relativeWeight;

    private double[] baseStationaryDistribution;
    private double[] storedBaseStationaryDistribution;

    private boolean stationaryDistributionKnown;
    private boolean storedStationaryDistributionKnown;
}
