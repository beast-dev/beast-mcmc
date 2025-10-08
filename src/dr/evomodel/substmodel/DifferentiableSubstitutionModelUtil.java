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

import java.util.Arrays;

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

        int[] complexIndices = new int[stateCount];
        int[] realIndices = new int[stateCount];
        int numComplexPairs = getComplexEigenValueFirstIndices(eigenValues, complexIndices);
        int numRealEigenValues = getRealEigenValueIndices(eigenValues, realIndices);

        getTripleMatrixMultiplication(stateCount, inverseEigenVectors, differentialMassMatrix, eigenVectors);

        setZeros(differentialMassMatrix);

        for (int i = 0; i < numRealEigenValues; i++) {
            final double eigenValueI = eigenValues[realIndices[i]];
            for (int j = 0; j < numRealEigenValues; j++) {
                final double eigenValueJ = eigenValues[realIndices[j]];
                if (i == j || Math.abs(eigenValueI - eigenValueJ) < threshold) {
                    differentialMassMatrix.set(realIndices[i], realIndices[j], differentialMassMatrix.get(realIndices[i], realIndices[j]) * time);
                } else {
                    differentialMassMatrix.set(realIndices[i], realIndices[j], differentialMassMatrix.get(realIndices[i], realIndices[j]) == 0 ? 0 : differentialMassMatrix.get(realIndices[i], realIndices[j]) * (1.0 - Math.exp((eigenValueJ - eigenValueI) * time)) / (eigenValueI - eigenValueJ));
                }
            }
        }

        for (int i = 0; i < numRealEigenValues; i++) {
            final double eigenValueI = eigenValues[realIndices[i]];
            final int iIndex = realIndices[i];
            for (int j = 0; j < numComplexPairs; j++) {
                final int jIndex = complexIndices[j];
                final double realEigenValue = eigenValues[jIndex];
                final double imagEigenValue = eigenValues[jIndex + stateCount];
                final double Vij = differentialMassMatrix.get(iIndex, jIndex);
                final double Vijp1 = differentialMassMatrix.get(iIndex, jIndex + 1);
                final double expSineIntegral = getExpSineIntegral(time, realEigenValue - eigenValueI, imagEigenValue);
                final double expCosineIntegral = getExpCosineIntegral(time, realEigenValue - eigenValueI, imagEigenValue);

                final double outij = Vij * expCosineIntegral - Vijp1 * expSineIntegral;
                final double outijp1 = Vij * expSineIntegral + Vijp1 * expCosineIntegral;

                differentialMassMatrix.set(iIndex, jIndex, outij);
                differentialMassMatrix.set(iIndex, jIndex + 1, outijp1);
            }
        }

        for (int i = 0; i < numComplexPairs; i++) {
            final int iIndex = complexIndices[i];
            final double realEigenValueI = eigenValues[iIndex];
            final double imagEigenValueI = eigenValues[iIndex + stateCount];
            final double cosineBt = Math.cos(imagEigenValueI * time);
            final double sineBt = Math.sin(imagEigenValueI * time);
            for (int j = 0; j < numRealEigenValues; j++) {
                final int jIndex = realIndices[j];
                final double realEigenValueJ = eigenValues[jIndex];

                final double Vij = differentialMassMatrix.get(iIndex, jIndex);
                final double Vip1j = differentialMassMatrix.get(iIndex + 1, jIndex);

                final double expCosineConvolution = getExpCosineConvolution(time, realEigenValueJ - realEigenValueI, imagEigenValueI, imagEigenValueI * time);
                final double expSineConvolution = getExpSineConvolution(time, realEigenValueJ - realEigenValueI, imagEigenValueI, imagEigenValueI * time);

                final double tmpIj = Vij * expCosineConvolution + Vip1j * expSineConvolution;
                final double tmpIp1j = - Vij * expSineConvolution + Vip1j * expCosineConvolution;

                differentialMassMatrix.set(iIndex, jIndex, cosineBt * tmpIj - sineBt * tmpIp1j);
                differentialMassMatrix.set(iIndex + 1, jIndex, sineBt * tmpIj + cosineBt * tmpIp1j);
            }

            for (int j = 0; j < numComplexPairs; j++) {
                final int jIndex = complexIndices[j];
                final double realEigenValueJ = eigenValues[jIndex];
                final double imagEigenValueJ = eigenValues[jIndex + stateCount];

                final double Vij = differentialMassMatrix.get(iIndex, jIndex);
                final double Vijp1 = differentialMassMatrix.get(iIndex, jIndex + 1);
                final double Vip1j = differentialMassMatrix.get(iIndex + 1, jIndex);
                final double Vip1jp1 = differentialMassMatrix.get(iIndex + 1, jIndex + 1);

                final double expCosineXPlusY = i == j ? time * Math.cos(imagEigenValueI * time) : getExpCosineConvolution(time, realEigenValueJ - realEigenValueI, imagEigenValueI - imagEigenValueJ, imagEigenValueI * time);
                final double expCosineXMinusY = i == j ? Math.sin(imagEigenValueI * time) / imagEigenValueI : getExpCosineConvolution(time, realEigenValueJ - realEigenValueI, imagEigenValueI + imagEigenValueJ, imagEigenValueI * time);
                final double expSineXPlusY = i == j ? time * Math.sin(imagEigenValueI * time) : getExpSineConvolution(time, realEigenValueJ - realEigenValueI, imagEigenValueI - imagEigenValueJ, imagEigenValueI * time);
                final double expSineXMinusY = i == j ? 0 : getExpSineConvolution(time, realEigenValueJ - realEigenValueI, imagEigenValueI + imagEigenValueJ, imagEigenValueI * time);

                final double tmpIJ = (Vij + Vip1jp1) * expCosineXPlusY + (Vij - Vip1jp1) * expCosineXMinusY + (Vip1j - Vijp1) * expSineXPlusY + (Vip1j + Vijp1) * expSineXMinusY;
                final double tmpIJp1 = (Vijp1 - Vip1j) * expCosineXPlusY + (Vijp1 + Vip1j) * expCosineXMinusY + (Vip1jp1 + Vij) * expSineXPlusY + (Vip1jp1 - Vij) * expSineXMinusY;
                final double tmpIp1J = (Vip1j - Vijp1) * expCosineXPlusY + (Vip1j + Vijp1) * expCosineXMinusY - (Vip1jp1 + Vij) * expSineXPlusY + (Vip1jp1 - Vij) * expSineXMinusY;
                final double tmpIp1Jp1 = (Vij + Vip1jp1) * expCosineXPlusY - (Vij - Vip1jp1) * expCosineXMinusY + (Vip1j - Vijp1) * expSineXPlusY - (Vip1j + Vijp1) * expSineXMinusY;

                final double outIJ = 0.5 * (cosineBt * tmpIJ - sineBt * tmpIp1J);
                final double outIJp1 = 0.5 * (cosineBt * tmpIJp1 - sineBt * tmpIp1Jp1);
                final double outIp1J = 0.5 * (sineBt * tmpIJ + cosineBt * tmpIp1J);
                final double outIp1Jp1 = 0.5 * (sineBt * tmpIJp1 + cosineBt * tmpIp1Jp1);

                differentialMassMatrix.set(iIndex, jIndex, outIJ);
                differentialMassMatrix.set(iIndex + 1, jIndex, outIp1J);
                differentialMassMatrix.set(iIndex, jIndex + 1, outIJp1);
                differentialMassMatrix.set(iIndex + 1, jIndex + 1, outIp1Jp1);
            }
        }

        getTripleMatrixMultiplication(stateCount, eigenVectors, differentialMassMatrix, inverseEigenVectors);

        double[] outputArray = new double[stateCount * stateCount];

        for (int i = 0, length = stateCount * stateCount; i < length; ++i) {
            outputArray[i] = differentialMassMatrix.get(i);
        }

        return outputArray;
    }

    private static int getComplexEigenValueFirstIndices(double[] eigenValues, int[] indices) {
        final int stateCount = eigenValues.length / 2;
        Arrays.fill(indices, -1);
        int currentIndex = 0;
        for (int i = 0; i < stateCount; ++i) {
            final double imagEigenValue = eigenValues[i + stateCount];
            if (imagEigenValue != 0) {
                indices[currentIndex++] = i;
                assert(eigenValues[i + 1 + stateCount] == -imagEigenValue);
                i++;
            }
        }
        return currentIndex;
    }

    private static int getRealEigenValueIndices(double[] eigenValues, int[] indices) {
        final int stateCount = eigenValues.length / 2;
        Arrays.fill(indices, -1);
        int currentIndex = 0;
        for (int i = 0; i < stateCount; ++i) {
            final double imagEigenValue = eigenValues[i + stateCount];
            if (imagEigenValue == 0) {
                indices[currentIndex++] = i;
            }
        }
        return currentIndex;
    }

    private static double getExpCosineIntegral(double time, double expRate, double cosRate) {
        final double denominator = expRate * expRate + cosRate * cosRate;
        final double expProduct = Math.exp(expRate * time);
        final double numerator = cosRate * expProduct * Math.sin(cosRate * time) + expRate * expProduct * Math.cos(cosRate * time) - expRate;
        return numerator / denominator;
    }

    private static double getExpSineIntegral(double time, double expRate, double sinRate) {
        final double denominator = expRate * expRate + sinRate * sinRate;
        final double expProduct = Math.exp(expRate * time);
        final double numerator = expRate * expProduct * Math.sin(sinRate * time) - sinRate * expProduct * Math.cos(sinRate * time) + sinRate;
        return numerator / denominator;
    }

    private static double getExpCosineConvolution(double time, double expRate, double cosRate, double convolveConst) {
        return Math.cos(convolveConst) * getExpCosineIntegral(time, expRate, cosRate) - Math.sin(convolveConst) * getExpSineIntegral(time, expRate, cosRate);
    }

    private static double getExpSineConvolution(double time, double expRate, double sinRate, double convolveConst) {
        return Math.sin(convolveConst) * getExpCosineIntegral(time, expRate, sinRate) - Math.cos(convolveConst) * getExpSineIntegral(time, expRate, sinRate);
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

    protected static void setupQDerivative(BaseSubstitutionModel substitutionModel, double[] differentialRates,
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
