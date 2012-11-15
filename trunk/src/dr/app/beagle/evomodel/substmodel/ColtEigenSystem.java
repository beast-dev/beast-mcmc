/*
 * ColtEigenSystem.java
 *
 * Copyright (C) 2002-2012 Alexei Drummond, Andrew Rambaut & Marc A. Suchard
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

package dr.app.beagle.evomodel.substmodel;

import cern.colt.matrix.linalg.Property;
import cern.colt.matrix.linalg.Algebra;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.DoubleMatrix2D;
import dr.math.matrixAlgebra.RobustEigenDecomposition;
import dr.math.matrixAlgebra.RobustSingularValueDecomposition;

/**
 * @author Marc Suchard
 */
public class ColtEigenSystem implements EigenSystem {

    public EigenDecomposition decomposeMatrix(double[][] matrix) {

        final int stateCount = matrix.length;

        RobustEigenDecomposition eigenDecomp = new RobustEigenDecomposition(
                new DenseDoubleMatrix2D(matrix),maxIterations);

        DoubleMatrix2D eigenV = eigenDecomp.getV();
        DoubleMatrix2D eigenVInv;

        if (checkConditioning) {                                                
            RobustSingularValueDecomposition svd;
            try {
                svd = new RobustSingularValueDecomposition(eigenV, maxIterations);
            } catch (ArithmeticException ae) {
                System.err.println(ae.getMessage());
                return getEmptyDecomposition(stateCount);
            }
            if (svd.cond() > maxConditionNumber) {
                return getEmptyDecomposition(stateCount);
            }
        }

        try {
            eigenVInv = alegbra.inverse(eigenV);
        } catch (IllegalArgumentException e) {
            return getEmptyDecomposition(stateCount);
        }

        double[][] Evec = eigenV.toArray();
        double[][] Ievc = eigenVInv.toArray();
        double[] Eval = getAllEigenValues(eigenDecomp);

        if (checkConditioning) {        
            for (int i = 0; i < Eval.length; i++) {
                if (Double.isNaN(Eval[i])  ||
                    Double.isInfinite(Eval[i])) {
                    return getEmptyDecomposition(stateCount);
                } else if (Math.abs(Eval[i]) < 1e-10) {
                    Eval[i] = 0.0;
                }
            }
        }

        double[] flatEvec = new double[stateCount * stateCount];
        double[] flatIevc = new double[stateCount * stateCount];

        for (int i = 0; i < Evec.length; i++) {
            System.arraycopy(Evec[i], 0, flatEvec, i * stateCount, stateCount);
            System.arraycopy(Ievc[i], 0, flatIevc, i * stateCount, stateCount);
        }

        return new EigenDecomposition(flatEvec, flatIevc, Eval);
    }

    protected double[] getAllEigenValues(RobustEigenDecomposition decomposition) {
        return decomposition.getRealEigenvalues().toArray();
    }

    protected double[] getEmptyAllEigenValues(int dim) {
        return new double[dim];
    }

    protected EigenDecomposition getEmptyDecomposition(int dim) {
        return new EigenDecomposition(
                new double[dim * dim],
                new double[dim * dim],
                getEmptyAllEigenValues(dim)
        );
    }


    protected boolean checkConditioning = true;
    protected int maxConditionNumber = 1000;
    protected int maxIterations = 1000;

    private static final double minProb = Property.DEFAULT.tolerance();
    private static final Algebra alegbra = new Algebra(minProb);

}
