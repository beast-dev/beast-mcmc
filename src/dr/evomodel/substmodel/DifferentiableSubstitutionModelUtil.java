/*
 * DifferentiableSubstitutionModelUtil.java
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

import dr.math.matrixAlgebra.ReadableMatrix;
import dr.math.matrixAlgebra.WrappedMatrix;
import dr.evomodel.substmodel.DifferentialMassProvider.DifferentialWrapper.WrtParameter;
import dr.math.matrixAlgebra.WrappedVector;

/**
 * @author Marc A. Suchard
 * @author Xiang Ji
 * @author Andrew Holbrook
 */
public class DifferentiableSubstitutionModelUtil {

    static double[] getApproximateDifferentialMassMatrix(double time,
                                                         WrappedMatrix differentialMassMatrix) {

        final int dim = differentialMassMatrix.getDim();

        double[] outputArray = new double[dim];

        for (int i = 0; i < dim; ++i) {
            outputArray[i] = time * differentialMassMatrix.get(i);
        }
        return outputArray;
    }

    static double[] getExactDifferentialMassMatrix(double time,
                                                   WrappedMatrix differentialMassMatrix,
                                                   EigenDecomposition eigenDecomposition) {
        
        final int stateCount = differentialMassMatrix.getMajorDim();

        double[] eigenValues = eigenDecomposition.getEigenValues();
        WrappedMatrix eigenVectors = new WrappedMatrix.Raw(eigenDecomposition.getEigenVectors(), 0, stateCount, stateCount);
        WrappedMatrix inverseEigenVectors = new WrappedMatrix.Raw(eigenDecomposition.getInverseEigenVectors(), 0, stateCount, stateCount);

        getTripleMatrixMultiplication(stateCount, inverseEigenVectors, differentialMassMatrix, eigenVectors);

        setZeros(differentialMassMatrix);

        for (int i = 0; i < stateCount; i++) {
            for (int j = 0; j < stateCount; j++) {
                if (i == j || eigenValues[i] == eigenValues[j]) {
                    differentialMassMatrix.set(i, j, differentialMassMatrix.get(i, j) * time);
                } else {
                    differentialMassMatrix.set(i, j, differentialMassMatrix.get(i, j) == 0 ? 0 : differentialMassMatrix.get(i, j) * (1.0 - Math.exp((eigenValues[j] - eigenValues[i]) * time)) / (eigenValues[i] - eigenValues[j]));
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

    static double[] getAffineDifferentialMassMatrix(double time,
                                                    WrappedMatrix differentialMassMatrix,
                                                    EigenDecomposition ed) {

        double[] differentials = getApproximateDifferentialMassMatrix(time, differentialMassMatrix);

        final int stateCount = differentialMassMatrix.getMajorDim();
        assert (stateCount == differentialMassMatrix.getMinorDim()) ;

        double[] correction = new double[stateCount * stateCount];

        int index = findZeroEigenvalueIndex(ed.getEigenValues(), stateCount);

        double[] eigenVectors = ed.getEigenVectors();
        double[] inverseEigenVectors = ed.getInverseEigenVectors();

        double[] qQPlus = getQQPlus(eigenVectors, inverseEigenVectors, index, stateCount);
        double[] oneMinusQPlusQ = getOneMinusQPlusQ(qQPlus, stateCount);

        double[] tmp = new double[stateCount * stateCount];
        multiply(tmp, oneMinusQPlusQ, differentials, 1.0, stateCount);
        multiply(correction, tmp, qQPlus, 1.0, stateCount);

        for (int i = 0; i < differentials.length; ++i) {
            differentials[i] -= correction[i];
        }

        return differentials;
    }
    
    private static void multiply(double[] result, double[] left, double[] right, double scale, int stateCount) {

        for (int i = 0; i < stateCount; ++i) {
            for (int j = 0; j < stateCount; ++j) {
                double entry = 0.0;
                for (int k = 0; k < stateCount; ++k) {
                    entry += left[i * stateCount + k] * right[k * stateCount + j];
                }
                result[i * stateCount + j] = scale * entry;
            }
        }
    }

    private static void setZeros(WrappedMatrix matrix) {
        for (int i = 0; i < matrix.getMinorDim(); i++) {
            for (int j = 0; j < matrix.getMinorDim(); j++) {
                if (Math.abs(matrix.get(i, j)) < threshold) {
                    matrix.set(i, j, 0);
                }
            }
        }
    }

    static final double threshold = 1E-10;  //TODO: very bad magic threshold number

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
        setupQDerivative(substitutionModel, differentialRates, differentialFrequencies, differentialMassMatrix);
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

    private static void setupQDerivative(BaseSubstitutionModel substitutionModel, double[] differentialRates,
                                         double[] differentialFrequencies, double[][] differentialMassMatrix) {
        if (substitutionModel instanceof ComplexSubstitutionModel) {
            int i, j, k = 0;
            final int stateCount = differentialFrequencies.length;
            for (i = 0; i < stateCount; i++) {
                for (j = i + 1; j < stateCount; j++) {
                    final double thisRate = differentialRates[k++];
                    differentialMassMatrix[i][j] = thisRate * differentialFrequencies[j];
                }
            }
            for (j = 0; j < stateCount; j++) {
                for (i = j + 1; i < stateCount; i++) {
                    final double thisRate = differentialRates[k++];
                    differentialMassMatrix[i][j] = thisRate * differentialFrequencies[j];
                }
            }
        } else {
            substitutionModel.setupQMatrix(differentialRates, differentialFrequencies, differentialMassMatrix);
        }
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


    private static int findZeroEigenvalueIndex(double[] eigenvalues, int stateCount) {
        for (int i = 0; i < stateCount; ++i) {
            if (eigenvalues[i] == 0) {
                return i;
            }
        }
        return -1;
    }

    private static double[] getOneMinusQPlusQ(double[] qQPlus, int stateCount) {
        double[] result = new double[stateCount * stateCount];

        for (int i = 0; i < stateCount; ++i) {
            for (int j = 0; j < stateCount; ++j) {
                if (i == j) {
                    result[index12(i, j, stateCount)] = 1.0 - qQPlus[index12(i, j, stateCount)];
                } else {
                    result[index12(i, j, stateCount)] = - qQPlus[index12(i, j, stateCount)];
                }
            }
        }

        return result;
    }

    public static double[] getQQPlus(double[] eigenVectors, double[] inverseEigenVectors, int index, int stateCount) {

        double[] result = new double[stateCount * stateCount];

        for (int i = 0; i < stateCount; ++i) {
            for (int j = 0; j < stateCount; ++j) {
                double sum = 0.0;
                for (int k = 0; k < stateCount; ++k) {
                    if (k != index) {
                        sum += eigenVectors[i * stateCount + k] * inverseEigenVectors[k * stateCount + j];
                    }
                }
                result[i * stateCount + j] = sum;
            }
        }

//        double[] reduced = new double[stateCount];
//        for (int j = 0; j < stateCount; ++j) {
//            double sum = 0.0;
//            for (int k = 0; k < stateCount; ++k) {
//                if (k != index) {
//                    sum += eigenVectors[k] * inverseEigenVectors[k * stateCount + j];
//                }
//            }
//            reduced[j] = sum;
//        }
//        reduced[0] -=1;

        // TODO Determine the stateCount unique values and just return them

        return result;
    }

    public static double[] getQQPlus(final double[] eigenVectors,
                                     final double[] inverseEigenVectors,
                                     final double[] eigenValues,
                                     final int stateCount) {

        double[] result = new double[stateCount * stateCount];

        for (int i = 0; i < stateCount; ++i) {
            for (int j = 0; j < stateCount; ++j) {
                double sum = 0.0;
                for (int k = 0; k < stateCount; ++k) {
                    if (eigenValues[k] != 0.0) {
                        sum += eigenVectors[i * stateCount + k] * inverseEigenVectors[k * stateCount + j];
                    }
                }
                result[i * stateCount + j] = sum;
            }
        }

        double[] reduced = new double[stateCount];
        for (int j = 0; j < stateCount; ++j) {
            double sum = 0.0;
            for (int k = 0; k < stateCount; ++k) {
                if (eigenValues[k] != 0.0) {
                    sum += eigenVectors[k] * inverseEigenVectors[k * stateCount + j];
                }
            }
            reduced[j] = sum;
        }
        reduced[0] -=1;

        // TODO Determine the stateCount unique values and just return them

        return result;
    }

//    public static double[] getQPlusQ(double[] eigenVectors, double[] inverseEigenVectors, int index, int stateCount) {
//
//        double[] result = new double[stateCount * stateCount];
//
//        for (int i = 0; i < stateCount; ++i) {
//            for (int j = 0; j < stateCount; ++j) {
//                double sum = 0.0;
//                for (int k = 0; k < stateCount; ++k) {
//                    if (k != index) {
//                        sum += inverseEigenVectors[i * stateCount + k] * eigenVectors[k * stateCount + j];
//                    }
//                }
//                result[i * stateCount + j] = sum;
//            }
//        }
//
//        return result;
//    }

    private static int index12(int i, int j, int stateCount) {
        return i * stateCount + j;
    }

    private static int index21(int i, int j, int stateCount) {
        return j * stateCount + i;
    }
}
