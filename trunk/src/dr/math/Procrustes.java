package dr.math;

/**
 * Procrustination function based on procrustes.r
 *
 * @author Andrew Rambaut
 * @author Marc Suchard
 * @version $Id$
 */
public class Procrustes {
    public final static double[][] procrustinate(double[][] X, double[][] Xstar, boolean allowTranslation, boolean allowDilation) {

        int n = X.length;
        int m = X[0].length;

       if (Xstar.length != n) {
           throw new IllegalArgumentException("X and Xstar do not have same number of rows");
       }
       if (Xstar[0].length != m) {
           throw new IllegalArgumentException("X and Xstar do not have same number of columns");
       }

       double[][] Xnew = new double[m][n];

       if (allowTranslation) {
//           J <- diag(n) - 1/n * matrix(1, n, n)
//           diag(n) : diagonal matrix of size n x n
//           matrix(1, n, n) : an n x n matrix with all entries = 1
//           for n = 3, J = {{1, -2/3, -2/3}, {-2/3, 1, -2/3}, {-2/3, -2/3, 1}}
       }
       else {
//           J <- diag(n)
//           see above

       }

//       C <- t(Xstar) %*% J %*% X
//        t(XStar) : transpose of Xstar
//        X %*% Y : X (matrix-multiply) Y

//       svd.out <- svd(C)
//       R <- svd.out$v %*% t(svd.out$u)
//       s <- 1

       if (allowDilation) {
//           mat1 <- t(Xstar) %*% J %*% X %*% R
//           mat2 <- t(X) %*% J %*% X
//           s.numer <- 0
//           s.denom <- 0
//           for (i in 1:m) {
//               s.numer <- s.numer + mat1[i, i]
           // mat1[i,i] : entry (i,i) of matrix mat1, remember R uses Fortran numbering (starting with 1)
//               s.denom <- s.denom + mat2[i, i]
//           }
//           s <- s.numer/s.denom
       }

//       tt <- matrix(0, m, 1)
        // an m x 1 matrix of all 0s

       if (allowTranslation) {
//           tt <- 1/n * t(Xstar - s * X %*% R) %*% matrix(1, n, 1)
       }

//       X.new <- s * X %*% R + matrix(tt, nrow(X), ncol(X), byrow = TRUE)
//        s * X  : element-wise multiply of all entries in X by s
//        nrow(X) : return the number of rows in X
//        ncol(X) : return the number of cols in X
//       return(list(X.new = X.new, R = R, tt = tt, s = s))

       return Xnew;
    }
}
