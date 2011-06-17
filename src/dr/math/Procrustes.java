package dr.math;

import dr.xml.Reference;
import no.uib.cipr.matrix.SVD;
import org.apache.commons.math.linear.*;

/**
 * Procrustination function based on procrustes.r
 *
 * @author Andrew Rambaut
 * @author Marc Suchard
 * @version $Id$
 */
public class Procrustes {
    public final static RealMatrix procrustinate(RealMatrix X, RealMatrix Xstar, boolean allowTranslation, boolean allowDilation) {

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
//           diag(n) : diagonal matrix of size n x n
//           matrix(1, n, n) : an n x n matrix with all entries = 1
//           for n = 3, J = {{1, -2/3, -2/3}, {-2/3, 1, -2/3}, {-2/3, -2/3, 1}}

            for (int i = 0; i < n; i++) {
                J.setEntry(i, i, 1);

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
//        t(XStar) : transpose of Xstar
//        X %*% Y : X (matrix-multiply) Y

        RealMatrix C = X.multiply(J.multiply(Xstar.transpose()));

//       svd.out <- svd(C)
//       R <- svd.out$v %*% t(svd.out$u)

        SingularValueDecomposition SVD = new SingularValueDecompositionImpl(C);
        RealMatrix R = SVD.getV().multiply(SVD.getUT());

//       s <- 1
        double s = 1.0;

        if (allowDilation) {
//           mat1 <- t(Xstar) %*% J %*% X %*% R
            RealMatrix mat1 = R.multiply(X.multiply(J.multiply(Xstar.transpose())));

//           mat2 <- t(X) %*% J %*% X
            RealMatrix mat2 = X.multiply(J.multiply(X.transpose()));

//           s.numer <- 0
//           s.denom <- 0
            double numer = 0.0;
            double denom = 0.0;

//           for (i in 1:m) {
//               s.numer <- s.numer + mat1[i, i]
            // mat1[i,i] : entry (i,i) of matrix mat1, remember R uses Fortran numbering (starting with 1)
//               s.denom <- s.denom + mat2[i, i]
//           }
            for (int i = 0; i < m; i++) {
                numer = numer + mat1.getEntry(i, i);
                denom = denom + mat2.getEntry(i, i);
            }
//           s <- s.numer/s.denom
           s = numer / denom;
        }

//       tt <- matrix(0, m, 1)
        // an m x 1 matrix of all 0s
        RealMatrix tt = new Array2DRowRealMatrix(m, 1);

        if (allowTranslation) {
//           tt <- 1/n * t(Xstar - s * X %*% R) %*% matrix(1, n, 1)
            tt = Xstar.subtract(R.multiply(X.scalarMultiply(s))).scalarMultiply(- 1.0 / n);
        }

//       X.new <- s * X %*% R + matrix(tt, nrow(X), ncol(X), byrow = TRUE)
//        s * X  : element-wise multiply of all entries in X by s
//        nrow(X) : return the number of rows in X
//        ncol(X) : return the number of cols in X

        // AR - I don't think this is correct (I am not sure what the byrow = TRUE does above):
        RealMatrix Xnew = R.multiply(X).scalarMultiply(s).add(tt);

//       return(list(X.new = X.new, R = R, tt = tt, s = s))

        return Xnew;
    }
}
