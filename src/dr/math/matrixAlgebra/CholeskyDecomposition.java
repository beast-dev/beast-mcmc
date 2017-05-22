/*
 * CholeskyDecomposition.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.math.matrixAlgebra;

/**
 * Created by IntelliJ IDEA.
 * User: msuchard
 * Date: Jan 12, 2007
 * Time: 9:05:44 PM
 * To change this template use File | Settings | File Templates.
 */
public class CholeskyDecomposition {

	/**
	 * Dimension of square matrix
	 */
	private int n;

	public boolean isSPD() {
		return isspd;
	}

	/**
	 * Symmetric and positive definite flag.
	 */
	private boolean isspd;

	public double[][] getL() {
		return L;
	}

	private double[][] L;

	public CholeskyDecomposition(Matrix A) throws IllegalDimension{
		this(A.components);
	}

	public CholeskyDecomposition(double[][] A) throws IllegalDimension {

		n = A.length;
		L = new double[n][n];
		isspd = (A[0].length == n);
		if (!isspd)
			throw new IllegalDimension("Cholesky decomposition is only defined for square matrices");
		// Main loop.
		for (int j = 0; j < n; j++) {
			double[] Lrowj = L[j];
			double d = 0.0;
			for (int k = 0; k < j; k++) {
				double[] Lrowk = L[k];
				double s = 0.0;
				for (int i = 0; i < k; i++) {
					s += Lrowk[i] * Lrowj[i];
				}
				Lrowj[k] = s = (A[j][k] - s) / L[k][k];
				d = d + s * s;
				isspd = isspd & (A[k][j] == A[j][k]);
			}
			d = A[j][j] - d;
			isspd = isspd & (d > 0.0);
			L[j][j] = Math.sqrt(Math.max(d, 0.0));
			/*for (int k = j+1; k < n; k++) {
						L[j][k] = 0.0;
					 }*/
		}

	}

    public static double[][] execute(double[] A, final int offset, final int n) {

        final double[][] L = new double[n][n];
//        boolean isspd = (A[0].length == n);
//        if (!isspd)
//            throw new IllegalDimension("Cholesky decomposition is only defined for square matrices");
        // Main loop.
        for (int j = 0; j < n; j++) {
            double[] Lrowj = L[j];
            double d = 0.0;
            for (int k = 0; k < j; k++) {
                double[] Lrowk = L[k];
                double s = 0.0;
                for (int i = 0; i < k; i++) {
                    s += Lrowk[i] * Lrowj[i];
                }
                Lrowj[k] = s = (A[offset + j * n + k] - s) / L[k][k];
                d = d + s * s;
//                isspd = isspd & (A[k][j] == A[j][k]);
            }
            d = A[offset + j * n + j] - d;
//            isspd = isspd & (d > 0.0);
            L[j][j] = Math.sqrt(Math.max(d, 0.0));
        }
        return L;
    }
}
