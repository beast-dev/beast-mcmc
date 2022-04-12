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

    static double[] getDifferentialMassMatrix(double time,
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

    private static void getTripleMatrixMultiplication(int stateCount, ReadableMatrix leftMatrix,
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

        final double[] differentialFrequencies = new double[stateCount];
        ((DifferentiableSubstitutionModel) substitutionModel).setupDifferentialFrequency(wrt, differentialFrequencies);

        double[][] differentialMassMatrix = new double[stateCount][stateCount];
        substitutionModel.setupQMatrix(differentialRates, differentialFrequencies, differentialMassMatrix);
        substitutionModel.makeValid(differentialMassMatrix, stateCount);

        final double weightedNormalizationGradient
                = ((DifferentiableSubstitutionModel) substitutionModel).getWeightedNormalizationGradient(
                        wrt, differentialMassMatrix, differentialFrequencies);

        for (int i = 0; i < stateCount; i++) {
            for (int j = 0; j < stateCount; j++) {
                differentialMassMatrix[i][j] -= Q[i * stateCount + j] * weightedNormalizationGradient;
            }
        }

        WrappedMatrix differential = new WrappedMatrix.ArrayOfArray(differentialMassMatrix);

        if (CHECK_COMMUTABILITY) {
            checkCommutability(differential, new WrappedMatrix.Raw(Q, 0, stateCount, stateCount));
        }

        return differential;
    }

    private static final boolean CHECK_COMMUTABILITY = false;
    private static final double COMMUTABILITY_CHECK_THRESHOLD = 0.01;

    public static boolean checkCommutability(WrappedMatrix x, WrappedMatrix y) {

        WrappedMatrix xy = product(x, y);
        WrappedMatrix yx = product(y, x);

//        System.err.println(xy);
//        System.err.println(yx);
//        System.err.println();

        boolean isCommutable = true;
        for (int i = 0; i < xy.getDim(); i++) {
            if (Math.abs(2.0 * (xy.get(i) - yx.get(i)) /  (xy.get(i) + yx.get(i)))> COMMUTABILITY_CHECK_THRESHOLD) {
                isCommutable = false;
            }
        }

//        if (isCommutable) {
//            System.err.println("Generator and its differential matrix commute.");
//        } else {
//            System.err.println("Generator and its differential matrix do not commute.");
//        }

        return isCommutable;
    }

    private static WrappedMatrix product(WrappedMatrix x, WrappedMatrix y) {
        final int majorDim = x.getMajorDim();
        final int minorDim = y.getMinorDim();

        final int innerDim = x.getMinorDim();

        if (innerDim != y.getMajorDim()) {
            return null;
        }

        WrappedMatrix result = new WrappedMatrix.Raw(new double[majorDim * minorDim], 0, majorDim, minorDim);

        for (int i = 0; i < majorDim; ++i) {
            for (int j = 0; j < minorDim; ++j) {
                double total = 0.0;
                for (int k = 0; k < innerDim; ++k) {
                    total += x.get(i, k) * y.get(k, j);
                }
                result.set(i, j, total);
            }
        }

        return result;
    }

}
