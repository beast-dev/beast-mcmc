/*
 * PolymorphismAwareSubstitutionModel.java
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

import dr.evolution.datatype.Codons;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.datatype.PolymorphismAwareDataType;
import dr.inference.model.*;
import dr.util.Citable;
import dr.util.Citation;

import java.util.List;

/**
 * PoMo model of sequence evolution
 *
 * @author Xiang Ji
 * @author Nicola De Maio
 * @author Ben Redelings
 * @author Marc A. Suchard
 */
public class PolymorphismAwareSubstitutionModel extends AbstractModel implements ActionEnabledSubstitution {

    private final SubstitutionModel baseSubstitutionModel;
    private final PolymorphismAwareDataType dataType;
    private final PolymorphismAwareFrequencyModel frequencyModel;
    private final int virtualPopSize;
    private final double[] baseInfinitesimalMatrix;
    private boolean baseSubstitutionKnown = false;

    public PolymorphismAwareSubstitutionModel(SubstitutionModel baseSubstitutionModel,
                                              PolymorphismAwareDataType dataType) {
        super("PoMoSubstitutionModel");
        this.baseSubstitutionModel = baseSubstitutionModel;

        this.frequencyModel = new PolymorphismAwareFrequencyModel(dataType, baseSubstitutionModel.getFrequencyModel().getFrequencyParameter(), baseSubstitutionModel);
        this.virtualPopSize = dataType.getVirtualPopSize();
        this.dataType = dataType;
        this.baseInfinitesimalMatrix = new double[dataType.getBaseDataType().getStateCount() * dataType.getBaseDataType().getStateCount()];
        syncBaseSubstitutionModel();
        addModel(baseSubstitutionModel);
    }


    @Override
    public int getNonZeroEntryCount() {
        int numberNucleotides;
        if (dataType.getBaseDataType() instanceof Nucleotides) {
            numberNucleotides = 1;
        } else if (dataType.getBaseDataType() instanceof Codons) {
            numberNucleotides = 3;
        } else {
            throw new RuntimeException("(PolymorphismAwareSubstitutionModel) Unsupported base data type: " + dataType.getBaseDataType().getClass().getName());
        }

        final int baseStateCount = dataType.getBaseDataType().getStateCount();
        final int originalNonZeros = baseStateCount * (3 * numberNucleotides + 1);
        final int polymorphismNonZeros = (baseStateCount * 3 * numberNucleotides + 1) / 2  * (virtualPopSize - 1) * 3;
        return originalNonZeros + polymorphismNonZeros;
    }

    @Override
    public EigenDecomposition getEigenDecomposition() {
        final int stateCount = dataType.getStateCount();

        if (eigenSystem == null) {
            eigenSystem = new ColtEigenSystem(stateCount);
        }

        double[][] infinitesimalMatrix = new double[stateCount][stateCount];
        int[] rowIndices = new int[getNonZeroEntryCount()];
        int[] colIndices = new int[getNonZeroEntryCount()];
        double[] values = new double[getNonZeroEntryCount()];
        getNonZeroEntries(rowIndices, colIndices, values);

        for (int i = 0; i < getNonZeroEntryCount(); i++) {
            infinitesimalMatrix[rowIndices[i]][colIndices[i]] = values[i];
        }

        return eigenSystem.decomposeMatrix(infinitesimalMatrix);
    }

    private EigenSystem eigenSystem = null;

    @Override
    public void getNonZeroEntries(int[] rowIndices, int[] colIndices, double[] values) {
        syncBaseSubstitutionModel();
        final int baseStateCount = dataType.getBaseDataType().getStateCount();
        double[] diagonal = new double[dataType.getStateCount()];
        int index = 0;
        for (int i = 0; i < baseStateCount; i++) {
            for (int j = 0; j < baseStateCount; j++) {
                final double baseRate = baseInfinitesimalMatrix[i * baseStateCount + j];
                if (i != j && baseRate != 0) {
                    rowIndices[index] = i;
                    final int colIndex = dataType.getState(new int[]{i, j}, new int[]{virtualPopSize - 1, 1});
                    colIndices[index] = colIndex;
                    values[index] = baseRate;
                    diagonal[i] -= baseRate;
                    index++;
                }
            }
        }

        for (int i = 0; i < baseStateCount - 1; i++) {
            for (int j = i + 1; j < baseStateCount; j++) {
                final double baseRate = baseInfinitesimalMatrix[i * baseStateCount + j];
                if (baseRate != 0) {
                    for (int k = 1; k < virtualPopSize; k++) {
                        final double rate = getDriftRate(k);
                        final int rowIndex = dataType.getState(new int[]{i, j}, new int[]{k, virtualPopSize - k});
                        final int plusOneColIndex = dataType.getState(new int[]{i, j}, new int[]{k + 1, virtualPopSize - k - 1});
                        final int minusOneColIndex = dataType.getState(new int[]{i, j}, new int[]{k - 1, virtualPopSize - k + 1});
                        rowIndices[index] = rowIndex;
                        colIndices[index] = plusOneColIndex;
                        values[index] = rate;
                        index++;
                        rowIndices[index] = rowIndex;
                        colIndices[index] = minusOneColIndex;
                        values[index] = rate;
                        index++;
                        diagonal[rowIndex] -= 2 * rate;
                    }
                }
            }
        }

        for (int i = 0; i < dataType.getStateCount(); i++) {
            if (diagonal[i] != 0) {
                rowIndices[index] = i;
                colIndices[index] = i;
                values[index] = diagonal[i];
                index++;
            }
        }
    }

    private double getDriftRate(int k) {
        return ((double) k * (virtualPopSize - k)) / (double) virtualPopSize;
    }

    private void syncBaseSubstitutionModel() {
        if (!baseSubstitutionKnown) {
            baseSubstitutionModel.getInfinitesimalMatrix(baseInfinitesimalMatrix);
            baseSubstitutionKnown = true;
        }
    }

    @Override
    public void getTransitionProbabilities(double distance, double[] matrix) {
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


    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == baseSubstitutionModel) {
            baseSubstitutionKnown = false;
            fireModelChanged();
        }
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {

    }


    @Override
    protected void storeState() {

    }

    @Override
    protected void restoreState() {
        baseSubstitutionKnown = false;
    }

    @Override
    protected void acceptState() {

    }

}
