package dr.evomodel.treedatalikelihood.continuous;

import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;

//public class MatrixExponential {
//}

public class MatrixExponential {

    private static final Algebra ALG = Algebra.DEFAULT;

    // Higham's theta_13 bound for 1-norm
    private static final double THETA_13 = 4.25;

    private MatrixExponential() {}

    /** expm(A) via scaling-and-squaring Padé(13). */
    public static DoubleMatrix2D expmPade13(DoubleMatrix2D A) {
        final int n = A.rows();
        if (n != A.columns()) {
            throw new IllegalArgumentException("A must be square");
        }
        // Quick exits
        if (n == 0) return new DenseDoubleMatrix2D(0, 0);
        if (isZero(A)) return eye(n);

        // 1-norm for scaling decision
        double A1 = norm1(A);

        int s = 0;
        if (A1 > THETA_13) {
            s = (int)Math.max(0, Math.ceil(log2(A1 / THETA_13)));
        }

        // A_s = A / 2^s
        DoubleMatrix2D As = A.copy();
        if (s > 0) scaleInPlace(As, 1.0 / Math.pow(2.0, s));

        // Powers: A2, A4, A6
        DoubleMatrix2D A2 = ALG.mult(As, As);
        DoubleMatrix2D A4 = ALG.mult(A2, A2);
        DoubleMatrix2D A6 = ALG.mult(A2, A4);

        // Higham's Padé(13) coefficients c0..c13 (c[0] is c0)
        final double[] c = new double[] {
                64764752532480000.0,
                32382376266240000.0,
                7771770303897600.0,
                1187353796428800.0,
                129060195264000.0,
                10559470521600.0,
                670442572800.0,
                33522128640.0,
                1323241920.0,
                40840800.0,
                960960.0,
                16380.0,
                182.0,
                1.0
        };

// Cache identity once
        DoubleMatrix2D I = eye(n);

// ---- Correct U and V per Higham ----

// Build tmp1 for the "odd" chain (c13, c11, c9)
        DoubleMatrix2D tmp1 = linComb(A6, c[13], A4, c[11]); // c13*A6 + c11*A4
        axpyInPlace(tmp1, A2, c[9]);                         // + c9*A2

// U = A * ( A6*tmp1 + c7*A6 + c5*A4 + c3*A2 + c1*I )
        DoubleMatrix2D polyU = ALG.mult(A6, tmp1);  // A6*tmp1
        axpyInPlace(polyU, A6, c[7]);               // + c7*A6
        axpyInPlace(polyU, A4, c[5]);               // + c5*A4
        axpyInPlace(polyU, A2, c[3]);               // + c3*A2
        axpyInPlace(polyU, I,  c[1]);               // + c1*I
        DoubleMatrix2D U = ALG.mult(As, polyU);     // multiply by A (scaled)

// Build inner for the "even" chain (c12, c10, c8)
        DoubleMatrix2D innerV = linComb(A6, c[12], A4, c[10]); // c12*A6 + c10*A4
        axpyInPlace(innerV, A2, c[8]);                         // + c8*A2

// V = A6*innerV + c6*A6 + c4*A4 + c2*A2 + c0*I
        DoubleMatrix2D V = ALG.mult(A6, innerV);  // A6*(c12*A6 + c10*A4 + c8*A2)
        axpyInPlace(V, A6, c[6]);                 // + c6*A6
        axpyInPlace(V, A4, c[4]);                 // + c4*A4
        axpyInPlace(V, A2, c[2]);                 // + c2*A2
        axpyInPlace(V, I,  c[0]);                 // + c0*I

// Solve (V - U) X = (V + U)
        DoubleMatrix2D VmU = V.copy(); axpyInPlace(VmU, U, -1.0);
        DoubleMatrix2D VpU = V.copy(); axpyInPlace(VpU, U,  1.0);
        DoubleMatrix2D X   = ALG.solve(VmU, VpU);

// Squaring
        for (int k = 0; k < s; k++) X = ALG.mult(X, X);
        return X;
    }

    // -------- helpers --------

    private static DoubleMatrix2D eye(int n) {
        DenseDoubleMatrix2D I = new DenseDoubleMatrix2D(n, n);
        for (int i = 0; i < n; i++) I.setQuick(i, i, 1.0);
        return I;
    }

    private static boolean isZero(DoubleMatrix2D A) {
        final int r = A.rows(), c = A.columns();
        for (int i = 0; i < r; i++)
            for (int j = 0; j < c; j++)
                if (A.getQuick(i, j) != 0.0) return false;
        return true;
    }

    /** 1-norm: max column sum of |A|. */
    private static double norm1(DoubleMatrix2D A) {
        final int n = A.rows();
        double max = 0.0;
        for (int j = 0; j < n; j++) {
            double s = 0.0;
            for (int i = 0; i < n; i++) s += Math.abs(A.getQuick(i, j));
            if (s > max) max = s;
        }
        return max;
    }

    private static double log2(double x) { return Math.log(x) / Math.log(2.0); }

    private static void scaleInPlace(DoubleMatrix2D A, double alpha) {
        final int r = A.rows(), c = A.columns();
        for (int i = 0; i < r; i++)
            for (int j = 0; j < c; j++)
                A.setQuick(i, j, A.getQuick(i, j) * alpha);
    }

    /** Z = a*X + b*Y (new matrix). */
    private static DoubleMatrix2D linComb(DoubleMatrix2D X, double a,
                                          DoubleMatrix2D Y, double b) {
        DenseDoubleMatrix2D Z = new DenseDoubleMatrix2D(X.rows(), X.columns());
        for (int i = 0; i < X.rows(); i++) {
            for (int j = 0; j < X.columns(); j++) {
                Z.setQuick(i, j, a * X.getQuick(i, j) + b * Y.getQuick(i, j));
            }
        }
        return Z;
    }

    /** X <- X + alpha * Y */
    private static void axpyInPlace(DoubleMatrix2D X, DoubleMatrix2D Y, double alpha) {
        final int r = X.rows(), c = X.columns();
        for (int i = 0; i < r; i++)
            for (int j = 0; j < c; j++)
                X.setQuick(i, j, X.getQuick(i, j) + alpha * Y.getQuick(i, j));
    }
}
