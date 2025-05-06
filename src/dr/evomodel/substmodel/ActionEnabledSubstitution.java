/*
 * ActionEnabledSubstitution.java
 *
 * Copyright (c) 2002-2018 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

import java.util.Arrays;

/**
 * @author Marc A. Suchard
 * @author Xiang Ji
 */
public interface ActionEnabledSubstitution extends SubstitutionModel {
    public int getNonZeroEntryCount();

    public void getNonZeroEntries(int[] rowIndices, int[] colIndices, double[] values);

    public class ActionEnabledSubstitutionWrap extends AbstractModel implements ActionEnabledSubstitution {

        private SubstitutionModel substitutionModel;

        private boolean substitutionModelKnown;

        private int[] rowIndices;

        private int[] colIndices;

        private double[] values;

        private double[] Q;

        private int stateCount;


        private int numNonZeroEntries;

        public ActionEnabledSubstitutionWrap(String name, SubstitutionModel substitutionModel) {
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
            substitutionModel.getInfinitesimalMatrix(matrix);
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
}


