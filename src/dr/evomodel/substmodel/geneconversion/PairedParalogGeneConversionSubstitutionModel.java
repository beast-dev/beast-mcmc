/*
 * PairedParalogGeneConversionSubstitutionModel.java
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

import dr.evolution.datatype.DataType;
import dr.evolution.datatype.PairedDataType;
import dr.evomodel.substmodel.*;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * @author Xiang Ji
 * @author Jeff Thorne
 * @author Marc A. Suchard
 */
public class PairedParalogGeneConversionSubstitutionModel extends AbstractModel implements Citable, ActionEnabledSubstitution {

    private final ActionEnabledSubstitution baseSubstitutionModel;
    private final Parameter igcRateParameter;
    protected PairedDataType dataType;

    private final RateCase rateCase;
    private final NumberParalog numberParalog;

    private final PairedParalogFrequencyModel frequencyModel;

    private final double[] values;

    private final int[] rowIndices;
    private final int[] colIndices;

    private boolean substitutionKnown;
    private int numNonZeros;

    public PairedParalogGeneConversionSubstitutionModel(String name,
                                                        BaseSubstitutionModel baseSubstitutionModel,
                                                        Parameter relativeGeneConversionRateParameter,
                                                        PairedDataType dataType,
                                                        RateCase rateCase,
                                                        NumberParalog numberParalog) {

        super(name);
        this.frequencyModel = new PairedParalogFrequencyModel(dataType, baseSubstitutionModel.getFrequencyModel().getFrequencyParameter());
        if (baseSubstitutionModel instanceof  ActionEnabledSubstitution) {
            this.baseSubstitutionModel = (ActionEnabledSubstitution) baseSubstitutionModel;
        } else {
            this.baseSubstitutionModel = new ActionEnabledSubstitutionWrap(baseSubstitutionModel.getModelName() + ".action.wrap", baseSubstitutionModel);
        }
        this.igcRateParameter = relativeGeneConversionRateParameter;
        this.dataType = dataType;
        this.rateCase = rateCase;
        this.numberParalog = numberParalog;

        this.values = new double[numberParalog.getMaxNonZeros(baseSubstitutionModel, frequencyModel, dataType)];
        this.rowIndices = new int[numberParalog.getMaxNonZeros(baseSubstitutionModel, frequencyModel, dataType)];
        this.colIndices = new int[numberParalog.getMaxNonZeros(baseSubstitutionModel, frequencyModel, dataType)];
        this.substitutionKnown = false;

        addVariable(relativeGeneConversionRateParameter);
        addModel(this.baseSubstitutionModel);
    }


    @Override
    public Citation.Category getCategory() {
        return Citation.Category.SUBSTITUTION_MODELS;
    }

    @Override
    public String getDescription() {
        return "Using pairwise interlocus gene conversion substitution model.";
    }

    @Override
    public List<Citation> getCitations() {
        return Collections.singletonList(CITATION);
    }

    private static final Citation CITATION = new Citation(
            new Author[]{
                    new Author("X", "Ji"),
                    new Author( "A", "Griffing"),
                    new Author("J", "Thorne"),
            },
            "A phylogenetic approach finds abundant interlocus gene conversion in yeast",
            2016,
            "Molecular Biology and Evolution",
            33,
            2469, 2476,
            Citation.Status.PUBLISHED);

    @Override
    public void getTransitionProbabilities(double distance, double[] matrix) {
        throw new RuntimeException("Not yet implemented!");
    }

    @Override
    public EigenDecomposition getEigenDecomposition() {
        throw new RuntimeException("Not yet implemented!");
    }

    @Override
    public FrequencyModel getFrequencyModel() {
        return frequencyModel;
    }

    @Override
    public void getInfinitesimalMatrix(double[] matrix) {
        throw new RuntimeException("Not yet implemented!");
    }

    @Override
    public DataType getDataType() {
        return dataType;
    }

    @Override
    public boolean canReturnComplexDiagonalization() {
        return false;
    }

    private void prepareSubstitutionModel() {
        Arrays.fill(rowIndices, 0);
        Arrays.fill(colIndices, 0);
        Arrays.fill(values, 0);

        numNonZeros = numberParalog.processSubstitutionModel(baseSubstitutionModel, igcRateParameter, dataType, rowIndices, colIndices, values);
        substitutionKnown = true;
    }

    @Override
    public int getNonZeroEntryCount() {
        if (!substitutionKnown) {
            prepareSubstitutionModel();
        }
        return numNonZeros;
    }

    @Override
    public void getNonZeroEntries(int[] outRowIndices, int[] outColIndices, double[] outValues) {
        if (!substitutionKnown) {
            prepareSubstitutionModel();
        }
        System.arraycopy(rowIndices, 0, outRowIndices, 0, numNonZeros);
        System.arraycopy(colIndices, 0, outColIndices, 0, numNonZeros);
        System.arraycopy(values, 0, outValues, 0, numNonZeros);
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        substitutionKnown = false;
        fireModelChanged();
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        substitutionKnown = false;
        fireModelChanged();
    }

    @Override
    protected void storeState() {

    }

    @Override
    protected void restoreState() {
        substitutionKnown = false;
    }

    @Override
    protected void acceptState() {

    }

    private class ActionEnabledSubstitutionWrap extends AbstractModel implements ActionEnabledSubstitution {

        private SubstitutionModel substitutionModel;

        private boolean substitutionModelKnown;

        private int[] rowIndices;

        private int[] colIndices;

        private double[] values;

        private double[] Q;

        private int stateCount;


        private int numNonZeroEntries;

        ActionEnabledSubstitutionWrap(String name, SubstitutionModel substitutionModel) {
            super(name);
            this.substitutionModel = substitutionModel;
            this.substitutionModelKnown = false;
            this.stateCount = substitutionModel.getFrequencyModel().getFrequencyCount();

            this.rowIndices = new int[stateCount * stateCount];
            this.colIndices = new int[stateCount * stateCount];
            this.values = new double[stateCount * stateCount];
            this.Q = new double[stateCount * stateCount];
            processSubstitutionModel();
            addModel(substitutionModel);
        }

        private void processSubstitutionModel() {
            Arrays.fill(rowIndices, 0);
            Arrays.fill(colIndices, 0);
            Arrays.fill(values, 0);

            substitutionModel.getInfinitesimalMatrix(Q);

            numNonZeroEntries = 0;
            for (int row = 0; row < stateCount; row++) {
                for (int col = 0; col < stateCount; col++) {
                    final double value = Q[row * stateCount + col];
                    if (value != 0) {
                        rowIndices[numNonZeroEntries] = row;
                        colIndices[numNonZeroEntries] = col;
                        values[numNonZeroEntries] = value;
                        numNonZeroEntries++;
                    }
                }
            }

            substitutionModelKnown = true;

        }

        @Override
        public int getNonZeroEntryCount() {
            if (!substitutionModelKnown) {
                processSubstitutionModel();
            }
            return numNonZeroEntries;
        }

        @Override
        public void getNonZeroEntries(int[] inRowIndices, int[] inColIndices, double[] inValues) {
            if (!substitutionModelKnown) {
                processSubstitutionModel();
            }
            System.arraycopy(rowIndices, 0, inRowIndices, 0, numNonZeroEntries);
            System.arraycopy(colIndices, 0, inColIndices, 0, numNonZeroEntries);
            System.arraycopy(values, 0, inValues, 0, numNonZeroEntries);
        }

        @Override
        public void getTransitionProbabilities(double distance, double[] matrix) {
            throw new RuntimeException("Not yet implemented!");
        }

        @Override
        public EigenDecomposition getEigenDecomposition() {
            throw new RuntimeException("Not yet implemented!");
        }

        @Override
        public FrequencyModel getFrequencyModel() {
            return substitutionModel.getFrequencyModel();
        }

        @Override
        public void getInfinitesimalMatrix(double[] matrix) {
            throw new RuntimeException("Not yet implemented!");
        }

        @Override
        public DataType getDataType() {
            throw new RuntimeException("Not yet implemented!");
        }

        @Override
        public boolean canReturnComplexDiagonalization() {
            throw new RuntimeException("Not yet implemented!");
        }


        @Override
        protected void handleModelChangedEvent(Model model, Object object, int index) {
            substitutionModelKnown = false;
            fireModelChanged();
        }

        @Override
        protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
            substitutionModelKnown = false;
            fireModelChanged();
        }

        @Override
        protected void storeState() {

        }

        @Override
        protected void restoreState() {
            substitutionModelKnown = false;
        }

        @Override
        protected void acceptState() {

        }
    }

    public enum RateCase {
        SINGLE("single") {
            @Override
            double getIGCRate(Parameter igcRateParameter, int donorParalogIndex) {
                return igcRateParameter.getParameterValue(0);
            }
        },
        ASYMMETRIC("asymmetric") {
            @Override
            double getIGCRate(Parameter igcRateParameter, int donorParalogIndex) {
                return igcRateParameter.getParameterValue(donorParalogIndex);
            }
        };

        private final String name;

        RateCase(String name) {
            this.name = name;
        }

        abstract double getIGCRate(Parameter igcRateParameter, int donorParalogIndex);

        public static RateCase factory(String match) {
            for (RateCase type : RateCase.values()) {
                if (match.equalsIgnoreCase(type.name)) {
                    return type;
                }
            }
            return null;
        }
    }

    public enum NumberParalog {
        SINGLE("1") {
            @Override
            void setupQMatrix(BaseSubstitutionModel baseSubstitutionModel, PairedParalogFrequencyModel freqModel,
                              PairedDataType dataType, double[][] matrix, RateCase rateCase, Parameter igcRateParameter) {

                final int baseNumStates = freqModel.getFrequencyParameter().getDimension();
                double[] infinitesimalMatrix = new double[baseNumStates * baseNumStates];
                baseSubstitutionModel.getInfinitesimalMatrix(infinitesimalMatrix);

                for (int i = 0; i < dataType.getStateCount(); i++) {
                    Arrays.fill(matrix[i], 0.0);
                }

                for (int i = 0; i < baseNumStates; i++) {

                    final int stateFrom = dataType.getState(i, i);

                    for (int j = 0; j < baseNumStates; j++) {
                        if (j != i) {
                            final int stateTo = dataType.getState(j, j);
                            matrix[stateFrom][stateTo] = infinitesimalMatrix[i * baseNumStates + j];
                        }
                    }
                }
            }

            @Override
            int processSubstitutionModel(ActionEnabledSubstitution baseSubstitution, Parameter igcRateParameter, PairedDataType dataType,
                                         int[] rowIndices, int[] colIndices, double[] values) {

                int[] baseRowIndices = new int[baseSubstitution.getNonZeroEntryCount()];
                int[] baseColIndices = new int[baseSubstitution.getNonZeroEntryCount()];

                baseSubstitution.getNonZeroEntries(baseRowIndices, baseColIndices, values);

                for (int i = 0; i < baseSubstitution.getNonZeroEntryCount(); i++) {
                    final int baseRowState = baseRowIndices[i];
                    final int baseColState = baseColIndices[i];

                    rowIndices[i] = dataType.getState(baseRowState, baseRowState);
                    colIndices[i] = dataType.getState(baseColState, baseColState);
                }

                return baseSubstitution.getNonZeroEntryCount();
            }

            @Override
            int getMaxNonZeros(BaseSubstitutionModel baseSubstitutionModel, PairedParalogFrequencyModel freqModel, PairedDataType dataType) {
                return baseSubstitutionModel.getFrequencyModel().getFrequencyCount() * baseSubstitutionModel.getFrequencyModel().getFrequencyCount();
            }
        },
        PAIR("2") {
            @Override
            void setupQMatrix(BaseSubstitutionModel baseSubstitutionModel, PairedParalogFrequencyModel freqModel,
                              PairedDataType dataType, double[][] matrix, RateCase rateCase, Parameter igcRateParameter) {
                final int baseNumStates = freqModel.getFrequencyParameter().getDimension();
                double[] infinitesimalMatrix = new double[baseNumStates * baseNumStates];
                baseSubstitutionModel.getInfinitesimalMatrix(infinitesimalMatrix);

                for (int i = 0; i < dataType.getStateCount(); i++) {
                    Arrays.fill(matrix[i], 0.0);
                    final int state1 = freqModel.getState1(i, baseNumStates);
                    final int state2 = freqModel.getState2(i, baseNumStates);

                    if (state1 != state2) {

                        for (int stateTo = 0; stateTo < baseNumStates; stateTo++) {

                            final int colIndex1 = dataType.getState(stateTo, state2);

                            if (stateTo != state1) {
                                matrix[i][colIndex1] = infinitesimalMatrix[state1 * baseNumStates + stateTo];

                                if (stateTo == state2) {
                                    matrix[i][colIndex1] += rateCase.getIGCRate(igcRateParameter, 1);
                                }
                            }

                            final int colIndex2 = dataType.getState(state1, stateTo);

                            if (stateTo != state2) {
                                matrix[i][colIndex2] = infinitesimalMatrix[state2 * baseNumStates + stateTo];

                                if (stateTo == state1) {
                                    matrix[i][colIndex2] += rateCase.getIGCRate(igcRateParameter, 0);
                                }
                            }
                        }
                    } else {
                        for (int stateTo = 0; stateTo < baseNumStates; stateTo++) {
                            if (stateTo != state1) {
                                final int colIndex1 = dataType.getState(stateTo, state2);
                                matrix[i][colIndex1] = infinitesimalMatrix[state1 * baseNumStates + stateTo];
                                final int colIndex2 = dataType.getState(state1, stateTo);
                                matrix[i][colIndex2] = infinitesimalMatrix[state2 * baseNumStates + stateTo];
                            }
                        }
                    }
                }
            }

            @Override
            int processSubstitutionModel(ActionEnabledSubstitution baseSubstitution, Parameter igcRateParameter, PairedDataType dataType,
                                         int[] rowIndices, int[] colIndices, double[] values) {

                final int stateCount = baseSubstitution.getFrequencyModel().getFrequencyCount();
                final double[] diagonal = new double[stateCount * stateCount];

                int[] baseRowIndices = new int[baseSubstitution.getNonZeroEntryCount()];
                int[] baseColIndices = new int[baseSubstitution.getNonZeroEntryCount()];
                double[] baseValues = new double[baseSubstitution.getNonZeroEntryCount()];

                baseSubstitution.getNonZeroEntries(baseRowIndices, baseColIndices, baseValues);


                int numNonZeros = 0;
                int numNonZeroDiagonal = 0;

                for (int i = 0; i < baseRowIndices.length; i++) {
                    final int fromState = baseRowIndices[i];
                    final int toState = baseColIndices[i];
                    final double rate = baseValues[i];

                    if (rate > 0) {
                        for (int j = 0; j < stateCount; j++) {
                            if (j != toState) {
                                rowIndices[numNonZeros] = dataType.getState(fromState, j);
                                colIndices[numNonZeros] = dataType.getState(toState, j);
                                values[numNonZeros] = rate;
                                diagonal[rowIndices[numNonZeros]] -= rate;
                                numNonZeros++;


                                rowIndices[numNonZeros] = dataType.getState(j, fromState);
                                colIndices[numNonZeros] = dataType.getState(j, toState);
                                values[numNonZeros] = rate;
                                diagonal[rowIndices[numNonZeros]] -= rate;
                                numNonZeros++;
                            }
                        }
                    } else {
                        numNonZeroDiagonal++;
                    }

                }

                final int numUniqueTransitions = baseRowIndices.length - numNonZeroDiagonal;

                final int numNoIGCTransitions = numUniqueTransitions * 2 * (stateCount - 1);
                final int numIGCTransitions = stateCount * (stateCount - 1) * 2;


                for (int i = 0; i < stateCount; i++) {
                    for (int j = 0; j < stateCount; j++) {
                        if (i != j) {
                            final int toState = dataType.getState(j, j);
                            rowIndices[numNonZeros] = dataType.getState(i, j);
                            colIndices[numNonZeros] = toState;
                            diagonal[rowIndices[numNonZeros]] -= igcRateParameter.getParameterValue(0);

                            numNonZeros++;

                            rowIndices[numNonZeros] = dataType.getState(j, i);
                            colIndices[numNonZeros] = toState;
                            diagonal[rowIndices[numNonZeros]] -= igcRateParameter.getParameterValue(0);
                            numNonZeros++;
                        }
                    }
                }

                Arrays.fill(values, numNoIGCTransitions, numNonZeros, igcRateParameter.getParameterValue(0));

                for (int i = 0; i < baseRowIndices.length; i++) {
                    final int fromState = baseRowIndices[i];
                    final int toState = baseColIndices[i];
                    final double rate = baseValues[i];

                    if (rate > 0) {
                        final int transitionIndex = fromState < toState ?
                                numNoIGCTransitions + (fromState * (stateCount - 1) + toState - 1) * 2 :
                                numNoIGCTransitions + (fromState * (stateCount - 1) + toState) * 2;

                        diagonal[rowIndices[transitionIndex]] -= rate;

                        assert(rowIndices[transitionIndex] == dataType.getState(fromState, toState) && colIndices[transitionIndex] == dataType.getState(toState, toState));

                        diagonal[rowIndices[transitionIndex + 1]] -= rate;

                        assert(rowIndices[transitionIndex + 1] == dataType.getState(toState, fromState) && colIndices[transitionIndex + 1] == dataType.getState(toState, toState));

                        values[transitionIndex] += rate;
                        values[transitionIndex + 1] += rate;
                    }
                }

                for (int i = 0; i < diagonal.length; i++) {
                    final double diagonalValue = diagonal[i];
                    if (diagonalValue != 0) {
                        rowIndices[numNonZeros] = i;
                        colIndices[numNonZeros] = i;
                        values[numNonZeros] = diagonalValue;
                        numNonZeros++;
                    }
                }

                return numNonZeros;
            }

            @Override
            int getMaxNonZeros(BaseSubstitutionModel baseSubstitutionModel, PairedParalogFrequencyModel freqModel, PairedDataType dataType) {
                final int totalPairedStateCount = freqModel.getFrequencyCount();
                final int baseStateCount = baseSubstitutionModel.getFrequencyModel().getFrequencyCount();
                return totalPairedStateCount * (2 * baseStateCount + 1);
            }
        };

        NumberParalog(String name) {
            this.name = name;
        }

        abstract void setupQMatrix(BaseSubstitutionModel baseSubstitutionModel, PairedParalogFrequencyModel freqModel,
                                   PairedDataType dataType, double[][] matrix, RateCase rateCase, Parameter igcRateParameter);

        abstract int processSubstitutionModel(ActionEnabledSubstitution baseSubstitution, Parameter igcRateParameter, PairedDataType dataType,
                                              int[] rowIndices, int[] colIndices, double[] values);
        private final String name;

        abstract int getMaxNonZeros(BaseSubstitutionModel baseSubstitutionModel, PairedParalogFrequencyModel freqModel,
                                    PairedDataType dataType);

        public static NumberParalog factory(String match) {
            for (NumberParalog type : NumberParalog.values()) {
                if (match.equalsIgnoreCase(type.name)) {
                    return type;
                }
            }
            return null;
        }
    }
}
