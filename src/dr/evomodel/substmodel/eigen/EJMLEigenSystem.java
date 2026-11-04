/*
 * EJMLEigenSystem.java
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
import org.ejml.alg.dense.decomposition.eig.SwitchingEigenDecomposition;
import org.ejml.data.DenseMatrix64F;
import static org.ejml.ops.CommonOps.invert;

/**
 * @author Xiang Ji
 * @author Marc Suchard
 */
public class EJMLEigenSystem implements EigenSystem {

    private final int stateCount;
    private final SwitchingEigenDecomposition ejmlEigenDecomposition;

    public EJMLEigenSystem(int stateCount) {
        this.stateCount = stateCount;
        this.ejmlEigenDecomposition = new SwitchingEigenDecomposition(stateCount * stateCount);

    }

    @Override
    public EigenDecomposition decomposeMatrix(double[][] matrix) {
        DenseMatrix64F A = new DenseMatrix64F(matrix);

        boolean decomposed = ejmlEigenDecomposition.decompose(A);

        int numberEigenValues = ejmlEigenDecomposition.getNumberOfEigenvalues();
        double[] eigenValues = new double[stateCount * 2];

        DenseMatrix64F eigenVector = new DenseMatrix64F(stateCount, stateCount);


        for (int i = 0; i < stateCount; i++) {
            eigenValues[i] = ejmlEigenDecomposition.getEigenvalue(i).getReal();
            if (ejmlEigenDecomposition.getEigenvalue(i).getImaginary() != 0) {
                eigenValues[i + stateCount] = ejmlEigenDecomposition.getEigenvalue(i).getImaginary();
            }

            for (int j = 0; j < stateCount; j++) {
                eigenVector.set(j, i, ejmlEigenDecomposition.getEigenVector(i).data[j]);
            }
        }

        double[] eigenVectors = eigenVector.getData().clone();

        invert(eigenVector);

        double[] inverseEigenVectors = eigenVector.getData().clone();

        return new EigenDecomposition(eigenVectors, inverseEigenVectors, eigenValues);
    }

    @Override
    public void computeExponential(EigenDecomposition ed, double time, double[] matrix) {
        throw new RuntimeException("Not yet implemented!");
    }

    @Override
    public double computeExponential(EigenDecomposition ed, double time, int i, int j) {
        throw new RuntimeException("Not yet implemented!");
    }
}
