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

import java.util.Arrays;

/**
 * @author Xiang Ji
 * @author Jeff Thorne
 * @author Marc A. Suchard
 */
public class SecondOrderMarkovSubstitutionModel extends BaseSubstitutionModel {
    
    private final BaseSubstitutionModel baseSubstitutionModel;
    private final SecondOrderMarkovPairedDataType pairedDataType;
    
    public SecondOrderMarkovSubstitutionModel(String name,
                                              SecondOrderMarkovPairedDataType dataType,
                                              BaseSubstitutionModel baseSubstitutionModel) {
        super(name, dataType, new SecondOrderMarkovFrequencyModel(SecondOrderMarkovFrequencyModel.NAME, baseSubstitutionModel, dataType));
        this.baseSubstitutionModel = baseSubstitutionModel;
        this.pairedDataType = (SecondOrderMarkovPairedDataType) getFrequencyModel().getDataType();
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
        super.setupMatrix();
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

                            matrix[fromState][toState] = infinitesimalMatrix[currentState * baseNumberStates + nextState];

                        }
                    }
                }
            }
        }

    }

    public static class SecondOrderMarkovFrequencyModel extends FrequencyModel {

        public static final String NAME = "SecondOrderMarkovFrequencyModel";
        private final BaseSubstitutionModel baseSubstitutionModel;

        private final SecondOrderMarkovPairedDataType dataType;

        private final int stateCount;

        private boolean frequencyKnown = false;
        private double[] frequencies;

        public SecondOrderMarkovFrequencyModel(String name,
                                               BaseSubstitutionModel substitutionModel,
                                               SecondOrderMarkovPairedDataType dataType) {
            super(name);
            this.baseSubstitutionModel = substitutionModel;
            this.dataType = dataType;
            this.stateCount = dataType.getStateCount();
            this.frequencies = new double[stateCount];
        }

        public double[] getFrequencies() {

            if (!frequencyKnown) {
                final double unif = 1.0 / stateCount;
                for (int i = 0; i < stateCount; i++) {
                    frequencies[i] = unif;
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
            if (baseDataType.isAmbiguousState(previousState) || baseDataType.isAmbiguousState(currentState)) {
                return getUnknownState();
            }
            return previousState * (baseDataType.getStateCount() - 1) + currentState;
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
