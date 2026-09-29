/*
 * Eigen3EigenSystem.java
 *
 * Copyright (c) 2002-2024 Alexei Drummond, Andrew Rambaut and Marc Suchard
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


package dr.evomodel.substmodel.eigen;

import dr.evomodel.substmodel.EigenDecomposition;
import dr.evomodel.substmodel.EigenSystem;

import java.util.Arrays;

/**
 * @author Xiang Ji
 * @author Marc Suchard
 */
public class Eigen3EigenSystem implements EigenSystem {

    private final int stateCount;
    private final EigenJNIWrapper eigenJNIWrapper;

    public Eigen3EigenSystem(int stateCount) {
        this.stateCount = stateCount;
        eigenJNIWrapper = new EigenJNIWrapper();
        eigenJNIWrapper.createInstance(1, stateCount);
    }

    @Override
    public EigenDecomposition decomposeMatrix(double[][] matrix) {

        final int stateCount = matrix.length;

        int[] indices = new int[2 * stateCount * stateCount];
        double[] values = new double[stateCount * stateCount];
        int nonZero = 0;

        for (int i = 0; i < stateCount; i++) {
            for (int j = 0; j < stateCount; j++) {
                final double value = matrix[i][j];
                if (value != 0) {
                    indices[2 * nonZero] = i;
                    indices[2 * nonZero + 1] = j;
                    values[nonZero] = value;
                    nonZero++;
                }
            }
        }
        eigenJNIWrapper.setMatrix(0, indices, values, nonZero);
        double[] eigenValues = new double[2 * stateCount];
        double[] eigenVectors = new double[stateCount * stateCount];
        double[] inverseEigenVectors = new double[stateCount * stateCount];
        eigenJNIWrapper.getEigenDecomposition(0, eigenValues, eigenVectors, inverseEigenVectors);

        return new EigenDecomposition(eigenVectors, inverseEigenVectors, eigenValues);
    }

    @Override
    public void computeExponential(EigenDecomposition eigen, double distance, double[] matrix) {
        double temp;

        if (eigen == null) {
            Arrays.fill(matrix, 0.0);
            return;
        }

        double[] Evec = eigen.getEigenVectors();
        double[] Eval = eigen.getEigenValues();
        double[] EvalImag = new double[stateCount];
        System.arraycopy(Eval, stateCount, EvalImag, 0, stateCount);
        double[] Ievc = eigen.getInverseEigenVectors();

        double[][] iexp = new double[stateCount][stateCount];

// Eigenvalues and eigenvectors of a real matrix A.
//
// If A is symmetric, then A = V*D*V' where the eigenvalue matrix D is diagonal
// and the eigenvector matrix V is orthogonal. I.e. A = V D V^t and V V^t equals
// the identity matrix.
//
// If A is not symmetric, then the eigenvalue matrix D is block diagonal with
// the real eigenvalues in 1-by-1 blocks and any complex eigenvalues,
// lambda + i*mu, in 2-by-2 blocks, [lambda, mu; -mu, lambda]. The columns
// of V represent the eigenvectors in the sense that A*V = V*D. The matrix
// V may be badly conditioned, or even singular, so the validity of the
// equation A = V D V^{-1} depends on the conditioning of V.

        for (int i = 0; i < stateCount; i++) {

            if (EvalImag[i] == 0) {
                // 1x1 block
                temp = Math.exp(distance * Eval[i]);
                for (int j = 0; j < stateCount; j++) {
                    iexp[i][j] = Ievc[i * stateCount + j] * temp;
                }
            } else {
                // 2x2 conjugate block
                // If A is 2x2 with complex conjugate pair eigenvalues a +/- bi, then
                // exp(At) = exp(at)*( cos(bt)I + \frac{sin(bt)}{b}(A - aI)).
                int i2 = i + 1;
                double b = EvalImag[i];
                double expat = Math.exp(distance * Eval[i]);
                double expatcosbt = expat * Math.cos(distance * b);
                double expatsinbt = expat * Math.sin(distance * b);

                for (int j = 0; j < stateCount; j++) {
                    iexp[i][j] = expatcosbt * Ievc[i * stateCount + j] +
                            expatsinbt * Ievc[i2 * stateCount + j];
                    iexp[i2][j] = expatcosbt * Ievc[i2 * stateCount + j] -
                            expatsinbt * Ievc[i * stateCount + j];
                }
                i++; // processed two conjugate rows
            }
        }

        int u = 0;
        for (int i = 0; i < stateCount; i++) {
            for (int j = 0; j < stateCount; j++) {
                temp = 0.0;
                for (int k = 0; k < stateCount; k++) {
                    temp += Evec[i * stateCount + k] * iexp[k][j];
                }
                matrix[u] = Math.abs(temp);
                u++;
            }
        }
    }

    @Override
    public double computeExponential(EigenDecomposition eigen, double distance, int i, int j) {
        throw new RuntimeException("Not yet implemented");
    }
}
