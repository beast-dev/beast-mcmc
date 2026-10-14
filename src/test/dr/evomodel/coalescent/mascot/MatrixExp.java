package test.dr.evomodel.coalescent.mascot;

/** Dense matrix exponential by scaling and squaring with a degree-13 Pade approximant (Higham 2005). */
public final class MatrixExp {

    private static final double[] B = {
            64764752532480000.0, 32382376266240000.0, 7771770303897600.0, 1187353796428800.0,
            129060195264000.0, 10559470521600.0, 670442572800.0, 33522128640.0, 1323241920.0,
            40840800.0, 960960.0, 16380.0, 182.0, 1.0};
    private static final double THETA13 = 5.371920351148152;

    private MatrixExp() {
    }

    public static double[][] expm(double[][] A0) {
        int n = A0.length;
        double norm = norm1(A0);
        int s = norm > THETA13 ? (int) Math.ceil(Math.log(norm / THETA13) / Math.log(2)) : 0;
        double scale = Math.pow(2, -s);
        double[][] A = scaled(A0, scale);
        double[][] I = identity(n);
        double[][] A2 = mul(A, A), A4 = mul(A2, A2), A6 = mul(A2, A4);
        double[][] U1 = lin(A6, B[13], A4, B[11], A2, B[9]);
        double[][] U2 = add(mul(A6, U1), lin(A6, B[7], A4, B[5], A2, B[3]));
        addInPlace(U2, I, B[1]);
        double[][] U = mul(A, U2);
        double[][] V1 = lin(A6, B[12], A4, B[10], A2, B[8]);
        double[][] V = add(mul(A6, V1), lin(A6, B[6], A4, B[4], A2, B[2]));
        addInPlace(V, I, B[0]);
        double[][] P = add(V, U);
        double[][] Q = sub(V, U);
        double[][] R = solve(Q, P);
        for (int k = 0; k < s; k++) R = mul(R, R);
        return R;
    }

    public static double[] rowTimes(double[] f, double[][] E) {
        int n = f.length;
        double[] out = new double[n];
        for (int x = 0; x < n; x++) for (int y = 0; y < n; y++) out[y] += f[x] * E[x][y];
        return out;
    }

    static double norm1(double[][] A) {
        int n = A.length;
        double best = 0;
        for (int j = 0; j < n; j++) {
            double s = 0;
            for (int i = 0; i < n; i++) s += Math.abs(A[i][j]);
            best = Math.max(best, s);
        }
        return best;
    }

    static double[][] identity(int n) {
        double[][] I = new double[n][n];
        for (int i = 0; i < n; i++) I[i][i] = 1;
        return I;
    }

    static double[][] scaled(double[][] A, double s) {
        int n = A.length;
        double[][] C = new double[n][n];
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) C[i][j] = A[i][j] * s;
        return C;
    }

    static double[][] mul(double[][] A, double[][] Bm) {
        int n = A.length;
        double[][] C = new double[n][n];
        for (int i = 0; i < n; i++) for (int k = 0; k < n; k++) {
            double a = A[i][k];
            if (a == 0) continue;
            for (int j = 0; j < n; j++) C[i][j] += a * Bm[k][j];
        }
        return C;
    }

    static double[][] add(double[][] A, double[][] Bm) {
        int n = A.length;
        double[][] C = new double[n][n];
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) C[i][j] = A[i][j] + Bm[i][j];
        return C;
    }

    static double[][] sub(double[][] A, double[][] Bm) {
        int n = A.length;
        double[][] C = new double[n][n];
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) C[i][j] = A[i][j] - Bm[i][j];
        return C;
    }

    static void addInPlace(double[][] A, double[][] Bm, double s) {
        int n = A.length;
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) A[i][j] += s * Bm[i][j];
    }

    static double[][] lin(double[][] X, double a, double[][] Y, double b, double[][] Z, double c) {
        int n = X.length;
        double[][] C = new double[n][n];
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) C[i][j] = a * X[i][j] + b * Y[i][j] + c * Z[i][j];
        return C;
    }

    /** Solve Q X = P by LU with partial pivoting. */
    static double[][] solve(double[][] Q0, double[][] P0) {
        int n = Q0.length;
        double[][] Q = new double[n][];
        double[][] P = new double[n][];
        for (int i = 0; i < n; i++) {
            Q[i] = Q0[i].clone();
            P[i] = P0[i].clone();
        }
        for (int col = 0; col < n; col++) {
            int piv = col;
            for (int r = col + 1; r < n; r++) if (Math.abs(Q[r][col]) > Math.abs(Q[piv][col])) piv = r;
            double[] t = Q[col];
            Q[col] = Q[piv];
            Q[piv] = t;
            t = P[col];
            P[col] = P[piv];
            P[piv] = t;
            double d = Q[col][col];
            for (int r = col + 1; r < n; r++) {
                double f = Q[r][col] / d;
                if (f == 0) continue;
                for (int c = col; c < n; c++) Q[r][c] -= f * Q[col][c];
                for (int c = 0; c < n; c++) P[r][c] -= f * P[col][c];
            }
        }
        double[][] X = new double[n][n];
        for (int r = n - 1; r >= 0; r--) {
            for (int c = 0; c < n; c++) {
                double s = P[r][c];
                for (int k = r + 1; k < n; k++) s -= Q[r][k] * X[k][c];
                X[r][c] = s / Q[r][r];
            }
        }
        return X;
    }
}
