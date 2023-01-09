/*
 * SecondOrderMarkovSubstitutionModel.java
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

package dr.evomodel.substmodel;

import dr.evolution.datatype.DataType;
import dr.evomodel.substmodel.spectra.SpectraJNIWrapper;
import dr.inference.model.Parameter;
import dr.xml.Reportable;

import java.util.Arrays;

/**
 * @author Xiang Ji
 * @author Jeff Thorne
 * @author Marc A. Suchard
 */
public class SecondOrderMarkovSubstitutionModel extends BaseSubstitutionModel implements Reportable {
    
    private final BaseSubstitutionModel baseSubstitutionModel;
    private final SecondOrderMarkovPairedDataType pairedDataType;

    private final ReversionRate reversionRate;

    private Boolean matrixKnown = false;
    
    public SecondOrderMarkovSubstitutionModel(String name,
                                              SecondOrderMarkovPairedDataType dataType,
                                              BaseSubstitutionModel baseSubstitutionModel,
                                              ReversionRate rate) {
        super(name, dataType, new SecondOrderMarkovFrequencyModel(SecondOrderMarkovFrequencyModel.NAME, dataType));
        this.baseSubstitutionModel = baseSubstitutionModel;
        this.pairedDataType = (SecondOrderMarkovPairedDataType) getFrequencyModel().getDataType();
        ((SecondOrderMarkovFrequencyModel) getFrequencyModel()).linkSecondOrderMarkovSubstitutionModel(this);
        this.reversionRate = rate;
    }

    protected void checkFrequencies() {

    }

    @Override
    protected void frequenciesChanged() {

    }

    @Override
    protected void ratesChanged() {

    }

    @Override
    protected void setupRelativeRates(double[] rates) {

    }

    @Override
    public boolean canReturnComplexDiagonalization() {
        return true;
    }

    @Override
    protected double setupMatrix() {
        if (!matrixKnown) {
            setupQMatrix(relativeRates, null, q);
        }
        return 1.0;
    }

    @Override
    protected EigenSystem getDefaultEigenSystem(int stateCount) {
        return new ComplexColtEigenSystem(stateCount);
    }

    @Override
    protected void setupQMatrix(double[] rates, double[] pi, double[][] matrix) {
        final int baseNumberStates = baseSubstitutionModel.stateCount;
        double[] infinitesimalMatrix = new double[baseNumberStates * baseNumberStates];
        baseSubstitutionModel.getInfinitesimalMatrix(infinitesimalMatrix);

        for (int previousState = 0; previousState < baseNumberStates; previousState++) {
            for (int currentState = 0; currentState < baseNumberStates; currentState++) {

                if (previousState != currentState) {
                    final int fromState = pairedDataType.getState(previousState, currentState);
                    Arrays.fill(matrix[fromState], 0);

                    for (int nextState = 0; nextState < baseNumberStates; nextState++) {
                        if (nextState != currentState) {
                            final int toState = pairedDataType.getState(currentState, nextState);

                            matrix[fromState][toState] = infinitesimalMatrix[currentState * baseNumberStates + nextState] + reversionRate.getReversionRate(fromState, toState);

                        }
                    }
                }
            }
        }
        matrixKnown = true;
    }

    public int getSparseQMatrix(int[] indices, double[] values, boolean transpose) {
        int nonZeroCount = 0;
        setupMatrix();
        for (int row = 0; row < getFrequencyModel().getFrequencyCount(); row++) {
            double rowSum = 0.0;
            for (int col = 0; col < getFrequencyModel().getFrequencyCount(); col++) {
                if (q[row][col] != 0) {
                    indices[nonZeroCount * 2] = transpose ? col : row;
                    indices[nonZeroCount * 2 + 1] = transpose ? row : col;
                    values[nonZeroCount] = q[row][col];
                    if (q[row][col] > 0)
                        rowSum += q[row][col];
                    nonZeroCount++;
                }
            }
            if (rowSum > 0) {
                indices[nonZeroCount * 2] = row;
                indices[nonZeroCount * 2 + 1] = row;
                values[nonZeroCount] = -rowSum;
                nonZeroCount++;
            }
        }
        return nonZeroCount;
    }

    public int getSparseEmbeddedMarkovChain(int[] indices, double[] values, boolean transpose) {
        int nonZeroCount = 0;
        setupMatrix();
        for (int row = 0; row < getFrequencyModel().getFrequencyCount(); row++) {
            double rowSum = 0.0;
            for (int col = 0; col < getFrequencyModel().getFrequencyCount(); col++) {
                if (row != col) rowSum += q[row][col];
            }
            for (int col = 0; col < getFrequencyModel().getFrequencyCount(); col++) {
                if (row != col && q[row][col] != 0) {
                    indices[nonZeroCount * 2] = transpose ? col : row;
                    indices[nonZeroCount * 2 + 1] = transpose ? row : col;
                    values[nonZeroCount] = q[row][col]/rowSum;
                    nonZeroCount++;
                }
            }
        }
        return nonZeroCount;
    }

    @Override
    public String getReport() {

        SpectraJNIWrapper spectra = SpectraJNIWrapper.loadLibrary();

        return spectra.getVersion();
    }

    public static class ReversionRate {
        private Parameter reversionRate;

        public ReversionRate(Parameter reversionRate) {
            this.reversionRate = reversionRate;
        }

        public double getReversionRate(int fromState, int toState) {
            return reversionRate.getParameterValue(0);
        }
    }

    public static class SecondOrderMarkovFrequencyModel extends FrequencyModel {

        public static final String NAME = "SecondOrderMarkovFrequencyModel";

        private SecondOrderMarkovSubstitutionModel secondOrderMarkovSubstitutionModel = null;

        private final SecondOrderMarkovPairedDataType dataType;

        private final int stateCount;

        private boolean frequencyKnown = false;
        private double[] frequencies;

        private final SpectraJNIWrapper spectra;

        public SecondOrderMarkovFrequencyModel(String name,
                                               SecondOrderMarkovPairedDataType dataType) {
            super(name);
            this.dataType = dataType;
            this.stateCount = dataType.getStateCount();
            this.frequencies = new double[stateCount];
            this.spectra = SpectraJNIWrapper.loadLibrary();
            spectra.createInstance(1, stateCount);
        }

        public void linkSecondOrderMarkovSubstitutionModel(SecondOrderMarkovSubstitutionModel substitutionModel) {
            this.secondOrderMarkovSubstitutionModel = substitutionModel;
        }

        public double[] getFrequencies() {

            if (!frequencyKnown) {

                if (secondOrderMarkovSubstitutionModel == null) {
                    throw new RuntimeException("Paired subsitution model unknown.");
                }

                double[] values = new double[dataType.getStateCount() * dataType.getStateCount()];
                int[] indices = new int[2 * dataType.getStateCount() * dataType.getStateCount()];

                int nonZeroCount = secondOrderMarkovSubstitutionModel.getSparseQMatrix(indices, values, true);
//                int nonZeroCount = secondOrderMarkovSubstitutionModel.getSparseEmbeddedMarkovChain(indices, values, false);

                spectra.setMatrix(0, indices, values, nonZeroCount);

                int numEigens = 1;
                double[] eigenValues = new double[2 * numEigens];
                double[] eigenVectors = new double[2 * stateCount * numEigens];
                spectra.getEigenVectors(0,numEigens , 1E-7, eigenValues, eigenVectors);

                double sum = 0;
                for (int i = 0; i < stateCount; i++) {
                    sum += eigenVectors[2 * i];
                }
                //TODO: check if all entries are of the same sign
                //TODO: check if eigen value is close to 0
                sum = sum == 0 ? 1.0 : sum;
                for (int i = 0; i < stateCount; i++) {
                    frequencies[i] = eigenVectors[2 * i] / sum;
                }
                frequencyKnown = true;
            }
            return frequencies;
        }

        public int getFrequencyCount() {
            return stateCount;
        }

        public void setFrequency(int i, double value) {
            // do nothing
        }

        public double getFrequency(int i) {
            if (! frequencyKnown) {
                getFrequencies();
            }
            return frequencies[i];
        }

        public DataType getDataType() {
            return dataType;
        }

    }

    public static class SecondOrderMarkovPairedDataType extends DataType {

        private final DataType baseDataType;

        public SecondOrderMarkovPairedDataType (DataType baseDataType) {
            this.baseDataType = baseDataType;
            stateCount = baseDataType.getStateCount() * (baseDataType.getStateCount() - 1);
            ambiguousStateCount = stateCount;
        }

        public final int getState(int previousState, int currentState) {
            if (baseDataType.isAmbiguousState(previousState) || baseDataType.isAmbiguousState(currentState) || previousState == currentState) {
                return getUnknownState();
            }
            return previousState < currentState ?
                    previousState * (baseDataType.getStateCount() - 1) + currentState - 1 :
                    previousState * (baseDataType.getStateCount() - 1) + currentState;
        }

        @Override
        public char[] getValidChars() {
            return null;
        }

        @Override
        public String getDescription() {
            return null;
        }

        @Override
        public int getType() {
            return -1;
        }
    }




}
