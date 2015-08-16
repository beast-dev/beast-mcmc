/*
 * Procrustes.java
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

package dr.math;

import org.apache.commons.math.linear.*;

/**
 * Procrustination function based on procrustes.r
 *
 * @author Andrew Rambaut
 * @author Marc Suchard
 * @version $Id$
 */
public class Procrustes {
    public Procrustes(RealMatrix X, RealMatrix Xstar, boolean allowTranslation, boolean allowDilation) {
        rowDimension = X.getRowDimension();
        columnDimension = X.getColumnDimension();

        if (Xstar.getRowDimension() != rowDimension) {
            throw new IllegalArgumentException("X and Xstar do not have the same number of rows");
        }
        if (Xstar.getColumnDimension() != columnDimension) {
            throw new IllegalArgumentException("X and Xstar do not have the same number of columns");
        }

        RealMatrix J = new Array2DRowRealMatrix(rowDimension, rowDimension);

        if (allowTranslation) {
//           J <- diag(n) - 1/n * matrix(1, n, n)
//           for n = 3, J = {{1, -2/3, -2/3}, {-2/3, 1, -2/3}, {-2/3, -2/3, 1}}

            for (int i = 0; i < rowDimension; i++) {
                J.setEntry(i, i, 1.0 - (1.0 / rowDimension));

                for (int j = i + 1; j < rowDimension; j++) {
                    J.setEntry(i, j, -1.0 / rowDimension);
                    J.setEntry(j, i, -1.0 / rowDimension);
                }
            }
        } else {
//           J <- diag(n)

            for (int i = 0; i < rowDimension; i++) {
                J.setEntry(i, i, 1);
            }

        }


//       C <- t(Xstar) %*% J %*% X

        RealMatrix C = Xstar.transpose().multiply(J.multiply(X));

//       svd.out <- svd(C)
//       R <- svd.out$v %*% t(svd.out$u)
//      NB: Apache math does a different SVD from R.  TODO Should use Colt library
        SingularValueDecomposition SVD = new SingularValueDecompositionImpl(C);
        R = SVD.getV().multiply(SVD.getUT());

//       s <- 1
        double s = 1.0; // scale = 1 unless dilation is being used

        if (allowDilation) {
//           mat1 <- t(Xstar) %*% J %*% X %*% R
            RealMatrix mat1 = Xstar.transpose().multiply(J.multiply(X.multiply(R)));

//           mat2 <- t(X) %*% J %*% X
            RealMatrix mat2 = X.transpose().multiply(J.multiply(X));

//           s.numer <- 0
//           s.denom <- 0
            double numer = 0.0;
            double denom = 0.0;

//           for (i in 1:m) {
//               s.numer <- s.numer + mat1[i, i]
//               s.denom <- s.denom + mat2[i, i]
//           }
            for (int i = 0; i < columnDimension; i++) {
                numer = numer + mat1.getEntry(i, i);
                denom = denom + mat2.getEntry(i, i);
            }
//           s <- s.numer/s.denom
            s = numer / denom;
        }
        this.s = s;

//       tt <- matrix(0, m, 1)
        RealMatrix tmpT = new Array2DRowRealMatrix(columnDimension, 1); // a translation vector of zero unless translation is being used

        if (allowTranslation) {
//           tt <- 1/n * t(Xstar - s * X %*% R) %*% matrix(1, n, 1)
            RealMatrix tmp = new Array2DRowRealMatrix(rowDimension, 1);
            for (int i = 0; i < rowDimension; i++) {
                tmp.setEntry(i, 0, 1);
            }
            tmpT = Xstar.subtract(X.multiply(R).scalarMultiply(s)).transpose().scalarMultiply(1.0 / rowDimension).multiply(tmp);
        }

        T = tmpT;
    }

    public final RealMatrix getTranslation() {
        return T.copy(); // NB Different from R
    }

    public final double getDilation() {
        return s;
    }

    public final RealMatrix getR() {
        return R.copy(); // NB Different from R
    }

    /**
     * procrustinate the complete matrix of coordinates
     * @param X the matrix containing coordinates (same dimensions as X in the constructor)
     * @return the transformed matrix
     */
    public final RealMatrix procrustinate(RealMatrix X) {
        if (X.getRowDimension() != rowDimension) {
            throw new IllegalArgumentException("X does not have the expected number of rows");
        }
        if (X.getColumnDimension() != columnDimension) {
            throw new IllegalArgumentException("X does not have the expected number of columns");
        }

//       X.new <- s * X %*% R + matrix(tt, nrow(X), ncol(X), byrow = TRUE)
        RealMatrix tt = new Array2DRowRealMatrix(rowDimension, columnDimension);

        for (int i = 0; i < rowDimension; i++) {
            tt.setRowMatrix(i, T.transpose());
        }

        // rotate, scale and translate
        return X.multiply(R).scalarMultiply(s).add(tt);  // Was a bug here
    }

    /**
     * procrustinate a single set of coordinates
     * @param X
     */
    public double[] procrustinate(double[] X) {
        if (X.length != columnDimension) {
            throw new IllegalArgumentException("X does not have the expected number of elements");
        }

        RealMatrix tmp = new Array2DRowRealMatrix(X);

        // rotate, scale and translate
        RealMatrix Xnew = tmp.multiply(R).scalarMultiply(s).add(T);

        return Xnew.getRow(0);
    }

    /**
     * procrustinate the complete matrix of coordinates
     */
    public final static RealMatrix procrustinate(RealMatrix X, RealMatrix Xstar, boolean allowTranslation, boolean allowDilation) {
       return new Procrustes(X, Xstar, allowTranslation, allowDilation).procrustinate(X);
    }


    private final RealMatrix R;
    private final RealMatrix T;
    private final double s;
    private final int rowDimension;
    private final int columnDimension;
}
