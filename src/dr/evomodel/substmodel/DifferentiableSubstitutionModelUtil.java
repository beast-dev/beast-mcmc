/*
 * DifferentiableSubstitutionModel.java
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

import dr.math.matrixAlgebra.ReadableMatrix;
import dr.math.matrixAlgebra.WrappedMatrix;
import dr.evomodel.substmodel.DifferentialMassProvider.DifferentialWrapper.WrtParameter;

/**
 * @author Marc A. Suchard
 * @author Xiang Ji
 */
public class DifferentiableSubstitutionModelUtil {

    public static double[] getDifferentialMassMatrix(double time,
                                                     int stateCount,
                                                     WrappedMatrix differentialMassMatrix,
                                                     EigenDecomposition eigenDecomposition) {

        double[] eigenValues = eigenDecomposition.getEigenValues();
        WrappedMatrix eigenVectors = new WrappedMatrix.Raw(eigenDecomposition.getEigenVectors(), 0, stateCount, stateCount);
        WrappedMatrix inverseEigenVectors = new WrappedMatrix.Raw(eigenDecomposition.getInverseEigenVectors(), 0, stateCount, stateCount);

        getTripleMatrixMultiplication(stateCount, inverseEigenVectors, differentialMassMatrix, eigenVectors);

        for (int i = 0; i < stateCount; i++) {
            for (int j = 0; j < stateCount; j++) {
                if (i == j || eigenValues[i] == eigenValues[j]) {
                    differentialMassMatrix.set(i, j, differentialMassMatrix.get(i, j) * time);
                } else {
                    differentialMassMatrix.set(i, j, differentialMassMatrix.get(i, j) * (1.0 - Math.exp((eigenValues[j] - eigenValues[i]) * time)) / (eigenValues[i] - eigenValues[j]));
                }
            }
        }

        getTripleMatrixMultiplication(stateCount, eigenVectors, differentialMassMatrix, inverseEigenVectors);

        double[] outputArray = new double[stateCount * stateCount];

        for (int i = 0, length = stateCount * stateCount; i < length; ++i) {
            outputArray[i] = differentialMassMatrix.get(i);
        }

        return outputArray;

    }

    public static void getTripleMatrixMultiplication(int stateCount, ReadableMatrix leftMatrix,
                                                      WrappedMatrix middleMatrix, ReadableMatrix rightMatrix) {

        double[][] tmpMatrix = new double[stateCount][stateCount];

        for (int i = 0; i < stateCount; i++) {
            for (int j = 0; j < stateCount; j++) {
                for (int k = 0; k < stateCount; k++) {
                    tmpMatrix[i][j] += middleMatrix.get(i, k) * rightMatrix.get(k, j);
                }
            }
        }

        for (int i = 0; i < stateCount; i++) {
            for (int j = 0; j < stateCount; j++) {
                double sumProduct = 0.0;
                for (int k = 0; k < stateCount; k++) {
                    sumProduct += leftMatrix.get(i, k) * tmpMatrix[k][j];
                }
                middleMatrix.set(i, j, sumProduct);
            }
        }
    }

    public static WrappedMatrix getInfinitesimalDifferentialMatrix(WrtParameter wrt,
                                                                   BaseSubstitutionModel substitutionModel) {
        if (!(substitutionModel instanceof DifferentiableSubstitutionModel)) {
            throw new RuntimeException("Not supported!");
        }
        final double normalizingConstant = substitutionModel.setupMatrix();

        final int stateCount = substitutionModel.getDataType().getStateCount();
        final int rateCount = substitutionModel.getRateCount(stateCount);

        final double[] Q = new double[stateCount * stateCount];
        substitutionModel.getInfinitesimalMatrix(Q);

        final double[] differentialRates = new double[rateCount];
        ((DifferentiableSubstitutionModel) substitutionModel).setupDifferentialRates(wrt, differentialRates, normalizingConstant);

        double[][] differentialMassMatrix = new double[stateCount][stateCount];
        substitutionModel.setupQMatrix(differentialRates, substitutionModel.getFrequencyModel().getFrequencies(), differentialMassMatrix);
        substitutionModel.makeValid(differentialMassMatrix, stateCount);

        final double weightedNormalizationGradient
                = ((DifferentiableSubstitutionModel) substitutionModel).getWeightedNormalizationGradient(differentialMassMatrix, substitutionModel.getFrequencyModel().getFrequencies());

        for (int i = 0; i < stateCount; i++) {
            for (int j = 0; j < stateCount; j++) {
                differentialMassMatrix[i][j] -= Q[i * stateCount + j] * weightedNormalizationGradient;
            }
        }

        return new WrappedMatrix.ArrayOfArray(differentialMassMatrix);
    }

}
