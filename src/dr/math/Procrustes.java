package dr.math;

import org.apache.commons.math.linear.*;

import java.awt.geom.AffineTransform;

/**
 * Procrustination function based on procrustes.r
 *
 * @author Andrew Rambaut
 * @author Marc Suchard
 * @version $Id$
 */
public class Procrustes {
    public Procrustes(RealMatrix X, RealMatrix Xstar, boolean allowTranslation, boolean allowDilation) {
        int n = X.getRowDimension();
        int m = X.getColumnDimension();

        if (Xstar.getRowDimension() != n) {
            throw new IllegalArgumentException("X and Xstar do not have same number of rows");
        }
        if (Xstar.getColumnDimension() != m) {
            throw new IllegalArgumentException("X and Xstar do not have same number of columns");
        }

        RealMatrix J = new Array2DRowRealMatrix(n, n);

        if (allowTranslation) {
//           J <- diag(n) - 1/n * matrix(1, n, n)
//           for n = 3, J = {{1, -2/3, -2/3}, {-2/3, 1, -2/3}, {-2/3, -2/3, 1}}

            for (int i = 0; i < n; i++) {
                J.setEntry(i, i, 1.0 - (1.0 / n));

                for (int j = i + 1; j < n; j++) {
                    J.setEntry(i, j, -1.0 / n);
                    J.setEntry(j, i, -1.0 / n);
                }
            }
        } else {
//           J <- diag(n)

            for (int i = 0; i < n; i++) {
                J.setEntry(i, i, 1);
            }

        }


//       C <- t(Xstar) %*% J %*% X

        RealMatrix C = Xstar.transpose().multiply(J.multiply(X));

//       svd.out <- svd(C)
//       R <- svd.out$v %*% t(svd.out$u)

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
            for (int i = 0; i < m; i++) {
                numer = numer + mat1.getEntry(i, i);
                denom = denom + mat2.getEntry(i, i);
            }
//           s <- s.numer/s.denom
            s = numer / denom;
        }
        this.s = s;

//       tt <- matrix(0, m, 1)
        RealMatrix tt = new Array2DRowRealMatrix(m, 1); // a translation vector of zero unless translation is being used

        if (allowTranslation) {
//           tt <- 1/n * t(Xstar - s * X %*% R) %*% matrix(1, n, 1)
            RealMatrix tmp = new Array2DRowRealMatrix(n, 1);
            for (int i = 0; i < n; i++) {
                tmp.setEntry(i, 0, 1);
            }
            tt = Xstar.subtract(X.multiply(R).scalarMultiply(s)).transpose().scalarMultiply(1.0 / n).multiply(tmp);
        }

//       X.new <- s * X %*% R + matrix(tt, nrow(X), ncol(X), byrow = TRUE)
        T = new Array2DRowRealMatrix(n, m);

        for (int i = 0; i < n; i++) {
            T.setRowMatrix(i, tt.transpose());
        }
    }

    /**
     * procrustinate the complete matrix of coordinates
     * @param X the matrix containing coordinates (same dimensions as X in the constructor)
     * @return the transformed matrix
     */
    public final RealMatrix procrustinate(RealMatrix X) {
        // rotate, scale and translate
        return X.multiply(R).scalarMultiply(s).add(T);
    }

    /**
     * procrustinate a single set of coordinates
     * @param X
     */
    public final void procrustinate(double[] X) {
        RealMatrix tmp = new Array2DRowRealMatrix(X);

        // rotate, scale and translate
        RealMatrix tmp2 = tmp.multiply(R).scalarMultiply(s).add(T);

        for (int i = 0; i < X.length; i++) {
            X[i] = tmp2.getEntry(i, 0);
        }
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
}
