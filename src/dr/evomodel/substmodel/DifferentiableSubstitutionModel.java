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

import dr.inference.model.Parameter;
import dr.math.matrixAlgebra.ReadableMatrix;
import dr.math.matrixAlgebra.WrappedMatrix;

/**
 * @author Marc A. Suchard
 * @author Xiang Ji
 */
public class DifferentiableSubstitutionModel{

    public static double[] getDifferentialMassMatrix(double time,
                                              Parameter parameter,
                                              int stateCount,
                                              WrappedMatrix.ArrayOfArray differentialMassMatrix,
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
        for (int i = 0; i < stateCount; i++) {
            System.arraycopy(differentialMassMatrix.getArrays()[i], 0, outputArray, i * stateCount, stateCount);
        }

        return outputArray;

    }

    public static void getTripleMatrixMultiplication(int stateCount, ReadableMatrix leftMatrix, WrappedMatrix middleMatrix, ReadableMatrix rightMatrix) {

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

}
