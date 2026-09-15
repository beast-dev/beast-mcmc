package dr.evomodel.coalescent;

public class GMRFDenseMatrixUtils {

    // computes aMat'*B
    public static double[][] transposeMultiply(double[][] aMat, double[][] bMat) {
        int n = aMat.length;
        int p1 = aMat[0].length;
        int p2 = bMat[0].length;
        double[][] result = new double[p1][p2];
        for (int i = 0; i < p1; i++) {
            for (int j = 0; j < p2; j++) {
                double sum = 0;
                for (int k = 0; k < n; k++) {
                    sum = sum + aMat[k][i]*bMat[k][j];
                }
                result[i][j] = sum;
            }
        }
        return result;
    }

    public static double[][] subtractMat(double[][] aMat, double[][] bMat) {
        int p = aMat.length;
        int q = aMat[0].length;
        double[][] result = new double[p][q];
        for (int i = 0; i < p; i++) {
            for (int j = 0; j < q; j++) {
                result[i][j] = aMat[i][j]-bMat[i][j];
            }
        }
        return result;
    }

    public static double[][] addMat(double[][] aMat, double[][] bMat) {
        int p = aMat.length;
        int q = aMat[0].length;
        double[][] result = new double[p][q];
        for (int i = 0; i < p; i++) {
            for (int j = 0; j < q; j++) {
                result[i][j] = aMat[i][j]+bMat[i][j];
            }
        }
        return result;
    }

    public static double[][] scaleMult(double[][] aMat, double c) {
        int p = aMat.length;
        int q = aMat[0].length;
        double[][] result = new double[p][q];
        for (int i = 0; i < p; i++) {
            for (int j = 0; j < q; j++) {
                result[i][j] = c*aMat[i][j];
            }
        }
        return result;
    }

    // Cholesky decomposition A = L*L' for symmetric positive definite
    // matrix aMat. Returns lower triangular factor L
    public static double[][] choleskyLower(double[][] aMat){
        int p = aMat.length;
        double[][] lMat = new double[p][p];
        for (int i = 0; i < p; i++) {
            for (int j = 0; j <= i; j++) {
                double sum = aMat[i][j];
                for (int k = 0; k < j; k++) {
                    sum = sum - lMat[i][k]*lMat[j][k];
                }
                if(i == j){
                    if(sum <= 0){
                        throw new IllegalArgumentException("Matrix is not positive definite." +
                                "Cholesky decomposition failed at at index " + i + ".");
                    }
                    lMat[i][j] = Math.sqrt(sum);
                }else{
                    lMat[i][j] = sum/lMat[j][j];
                }
            }
        }
        return lMat;
    }

    // Solves A*x = b for symmetric positive definite matrix A
    public static double[] solveSPD(double[][] aMat, double[] b){
        double[][] lMat = choleskyLower(aMat);
        int p = lMat.length;
        // lMat*y = b
        double[] y = new double[p];
        for (int i = 0; i < p; i++) {
            double sum = b[i];
            for (int k = 0; k < i; k++) {
                sum = sum - lMat[i][k]*y[k];
            }
            y[i] = sum/lMat[i][i];
        }
        // lMat'*x = y
        double[] x = new double[p];
        for (int i = p-1; i>=0; i--) {
            double sum = y[i];
            for (int k = i+1; k < p; k++) {
                sum = sum - lMat[k][i]*x[k];
            }
            x[i] = sum/lMat[i][i];
        }
        return x;
    }

    // computes log(det(aMat)) for symmetric positive definite aMat
    public static double logDetSPD(double[][] aMat) {
        double[][] lMat = choleskyLower(aMat);
        double logDet = 0;
        for (int i = 0; i < lMat.length; i++) {
            logDet = logDet + Math.log(lMat[i][i]);
        }
        return 2.0*logDet;
    }

    // Computes v'*A*v
    public static double quadraticForm(double[][] aMat, double[] v){
        int p = v.length;
        double result = 0;
        for (int i = 0; i < p; i++) {
            double sum = 0.0;
            for (int j = 0; j < p; j++) {
                sum = sum + aMat[i][j]*v[j];
            }
            result = result + v[i]*sum;
        }
        return result;
    }

    public static double[] diagonalOfInverseSPD(double[][] aMat){
        int p = aMat.length;
        double[] diag = new double[p];
        for (int i = 0; i < p; i++) {
            double[] e = new double[p];
            e[i] = 1.0;
            double[] x = solveSPD(aMat, e);
            diag[i] = x[i];
        }
        return diag;
    }

    public static double[] transposeMultiplyVector(double[][] aMat, double[] v){
        int n = aMat.length;
        int p = aMat[0].length;
        double[] result = new double[p];
        for (int j = 0; j < p; j++) {
            double sum = 0.0;
            for (int i = 0; i < n; i++) {
                sum = sum + aMat[i][j]*v[i];
            }
            result[j] = sum;
        }
        return result;
    }

    public static double[] multMatVec(double[][] aMat, double[] v){
        int rows = aMat.length;
        double[] result = new double[rows];
        for (int i = 0; i < rows; i++) {
            double sum = 0.0;
            for (int j = 0; j < v.length; j++) {
                sum = sum + aMat[i][j]*v[j];
            }
            result[i] = sum;
        }
        return result;
    }

    public static double[][] identityMat(int p){
        double[][] result = new double[p][p];
        for (int i = 0; i < p; i++) {
            result[i][i] = 1;
        }
        return result;
    }
}
