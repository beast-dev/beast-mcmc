/*
 * TwoParalogGeneConversionSubstitutionModel.java
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
import dr.evomodel.substmodel.BaseSubstitutionModel;
import dr.evomodel.substmodel.ComplexColtEigenSystem;
import dr.evomodel.substmodel.EigenSystem;
import dr.inference.model.Parameter;
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
public class PairedParalogGeneConversionSubstitutionModel extends BaseSubstitutionModel implements Citable {

    private final BaseSubstitutionModel baseSubstitutionModel;
    private final Parameter igcRateParameter;
    protected PairedDataType dataType;

    private final RateCase rateCase;
    private final NumberParalog numberParalog;

    public PairedParalogGeneConversionSubstitutionModel(String name,
                                                        BaseSubstitutionModel baseSubstitutionModel,
                                                        Parameter relativeGeneConversionRateParameter,
                                                        PairedDataType dataType,
                                                        RateCase rateCase,
                                                        NumberParalog numberParalog) {

        super(name, dataType, new PairedParalogFrequencyModel(dataType, baseSubstitutionModel.getFrequencyModel().getFrequencyParameter()));
        this.baseSubstitutionModel = baseSubstitutionModel;
        this.igcRateParameter = relativeGeneConversionRateParameter;
        this.dataType = dataType;
        this.rateCase = rateCase;
        this.numberParalog = numberParalog;

        addVariable(relativeGeneConversionRateParameter);
        addModel(baseSubstitutionModel);


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
        numberParalog.setupQMatrix(baseSubstitutionModel, (PairedParalogFrequencyModel) freqModel,
                dataType, matrix, rateCase, igcRateParameter);
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
        };

        NumberParalog(String name) {
            this.name = name;
        }

        abstract void setupQMatrix(BaseSubstitutionModel baseSubstitutionModel, PairedParalogFrequencyModel freqModel,
                                   PairedDataType dataType, double[][] matrix, RateCase rateCase, Parameter igcRateParameter);
        private final String name;

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
