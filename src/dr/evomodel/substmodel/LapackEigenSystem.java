/*
 * ColtEigenSystem.java
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

//import flanagan.complex.ComplexMatrix;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.EVD;

import java.util.Arrays;

/**
 * @author Marc Suchard
 */

public class LapackEigenSystem implements EigenSystem {

    public LapackEigenSystem(int stateCount) {
        this.stateCount = stateCount;
    }

    private static final boolean SORT = false;
    private static final boolean USE_LEFT = false;

    class SortedEvd {

        private final double[] realEigenvalues;
        private final double[] imagEigenvalues;

        private final DenseMatrix rightEigenvectors;
        private final DenseMatrix leftEigenvectors;

        private DenseMatrix inverseRightEigenvectors;

        SortedEvd(double[][] matrix) {

            EVD evd;
            try {
                evd = new EVD(matrix.length, USE_LEFT, true);
                evd.factor(new DenseMatrix(matrix));

            } catch (Exception exc) {
                throw new RuntimeException(exc.getMessage());
            }

            if (!SORT) {
                realEigenvalues = evd.getRealEigenvalues();
                imagEigenvalues = evd.getImaginaryEigenvalues();
                rightEigenvectors = evd.getRightEigenvectors();
                leftEigenvectors = evd.getLeftEigenvectors();
            } else {
                Integer[] indices = new Integer[stateCount];
                for (int i = 0; i < stateCount; ++i) {
                    indices[i] = i;
                }

                double[] rev = evd.getRealEigenvalues();
                double[] iev = evd.getImaginaryEigenvalues();

                Arrays.sort(indices, (lhs, rhs) -> -Double.compare(
                        rev[lhs] * rev[lhs] + iev[lhs] * iev[lhs],
                        rev[rhs] * rev[rhs] + iev[rhs] * iev[rhs]
                ));

                realEigenvalues = new double[stateCount];
                imagEigenvalues = new double[stateCount];
                rightEigenvectors = new DenseMatrix(stateCount, stateCount);
                leftEigenvectors = new DenseMatrix(stateCount, stateCount);

                for (int i = 0; i < stateCount; ++i) {
                    realEigenvalues[i] = rev[indices[i]];
                    imagEigenvalues[i] = iev[indices[i]];

                    for (int j = 0; j < stateCount; ++j) {
                        rightEigenvectors.set(i, j, evd.getRightEigenvectors().get(i, indices[j]));
                        if (USE_LEFT) {
                            leftEigenvectors.set(i, j, evd.getLeftEigenvectors().get(indices[i], j));
                        }
                    }
                }
            }
        }

        double[] getRealEigenvalues() { return realEigenvalues; }

        double[] getImaginaryEigenvalues() { return imagEigenvalues; }

        DenseMatrix getRightEigenvectors() { return rightEigenvectors; }

        @SuppressWarnings("unused")
        DenseMatrix getLeftEigenvectors() { return leftEigenvectors; }

        DenseMatrix getScaledInverseRightEigenvectors() {

            if (inverseRightEigenvectors == null) {
                DenseMatrix B = new DenseMatrix(stateCount, stateCount);
                for (int i = 0; i < stateCount; ++i) {
                        B.set(i, i, 1.0);
                }

                inverseRightEigenvectors = new DenseMatrix(stateCount, stateCount);
                rightEigenvectors.solve(B, inverseRightEigenvectors);
            }

            return inverseRightEigenvectors;
        }
    }

    public EigenDecomposition decomposeMatrix(double[][] matrix) {

        SortedEvd evd = new SortedEvd(matrix);

        // Flatten
        double[] compactEvec = new double[stateCount * stateCount];
        double[] compactIevc = new double[stateCount * stateCount];
        for (int i = 0; i < stateCount; ++i) {
            for (int j = 0; j < stateCount; ++j) {
                compactEvec[i * stateCount + j] = evd.getRightEigenvectors().get(i, j);
                compactIevc[i * stateCount + j] = evd.getScaledInverseRightEigenvectors().get(i, j);
            }
        }

        double[] compactEval = new double[stateCount * 2];
        System.arraycopy(evd.getRealEigenvalues(), 0, compactEval, 0, stateCount);
        System.arraycopy(evd.getImaginaryEigenvalues(), 0, compactEval, stateCount, stateCount);

        if (OLD_VALUES) {

            double[] lEval = new double[stateCount * 2];
            double[] lEvec = new double[stateCount * stateCount];
            double[] lIevc = new double[stateCount * stateCount];

            ComplexColtEigenSystem es = new ComplexColtEigenSystem(stateCount);
            EigenDecomposition old = es.decomposeMatrix(matrix);
            System.arraycopy(old.getEigenValues(), 0, lEval, 0, stateCount * 2);
            System.arraycopy(old.getEigenVectors(), 0, lEvec, 0, stateCount * stateCount);
            System.arraycopy(old.getInverseEigenVectors(), 0, lIevc, 0, stateCount * stateCount);

            return new LapackEigenDecomposition(lEvec, lIevc, lEval, matrix);

        } else {
            return new LapackEigenDecomposition(compactEvec, compactIevc, compactEval, matrix);
        }
    }

    private static final boolean OLD_VALUES = false;

    public double computeExponential(EigenDecomposition eigen, double distance, int i, int j) {
        throw new RuntimeException("Not yet implemented");
    }

    public void computeExponential(EigenDecomposition eigen, double distance, double[] matrix) {
        ComplexColtEigenSystem ces = new ComplexColtEigenSystem(stateCount);
        ces.computeExponential(eigen, distance, matrix);
    }

    @SuppressWarnings("comment")
/*
    public static void computeExponential(ComplexMatrix v, ComplexMatrix vInverse, org.hipparchus.complex.Complex[] eigenvalues,
                                          double distance, double[] matrix, int stateCount) {

        ComplexMatrix expD = new ComplexMatrix(stateCount, stateCount);

        for (int i = 0; i < stateCount; ++i) {
            org.hipparchus.complex.Complex expValue = eigenvalues[i].multiply(distance).exp();

            expD.setElement(i, i, expValue.getReal(), expValue.getImaginary());
        }

        ComplexMatrix result = v.times(expD.times(vInverse));

        for (int i = 0; i < stateCount; ++i) {
            for (int j = 0; j < stateCount; ++j) {
                matrix[i * stateCount + j] = result.getElementReference(i, j).getReal();
            }
        }
    }
*/

    protected final int stateCount;
}
