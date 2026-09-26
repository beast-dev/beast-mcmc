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

    // Solves LL'x = b given an already computed lower cholesky factor L
    public static double[] solveGivenCholeskyLower(double[][] lMat, double[] b){
        int p = lMat.length;
        double[] y = new double[p];
        for (int i = 0; i < p; i++) {
            double sum = b[i];
            for(int k = 0; k < i; k++){
                sum = sum - lMat[i][k]*y[k];
            }
            y[i] = sum/lMat[i][i];
        }
        double[] x = new double[p];
        for (int i = p-1; i>=0; i--) {
            double sum = y[i];
            for(int k = i+1; k < p; k++){
                sum = sum - lMat[k][i]*x[k];
            }
            x[i] = sum/lMat[i][i];
        }
        return x;
    }


    public static double logMvnDensityGivenPrecision(double[] x,
                                                     double[] mean,
                                                     double[][] precisionMat,
                                                     double logDetPrecision){
        int p = x.length;
        double[] diff = new double[p];
        for (int i = 0; i < p; i++) {
            diff[i] = x[i] - mean[i];
        }
        double quad = quadraticForm(precisionMat, diff);
        return -0.5*p*Math.log(2*Math.PI)+0.5*logDetPrecision-0.5*quad;
    }

    public static double[][] invertSPD(double[][] aMat){
        int p = aMat.length;
        double[][] result = new double[p][p];
        for (int j = 0; j < p; j++) {
            double[] e = new double[p];
            e[j] = 1.0;
            double[] col = solveSPD(aMat, e);
            for(int i=0; i<p; i++){
                result[i][j] = col[i];
            }
        }
        return result;
    }

    public static double logSumExp(double[] logValues){
        double max = Double.NEGATIVE_INFINITY;
        for(double v : logValues){
            if(v > max){
                max = v;
            }
        }
        if(Double.isInfinite(max)){
            return max;
        }
        double sum = 0.0;
        for (double v : logValues){
            sum += Math.exp(v-max);
        }
        return max + Math.log(sum);
    }

    // normalize the log weights into weights that sum to 1
    public static double[] normalizeLogWeights(double[] logWeights){
        double lse = logSumExp(logWeights);
        double[] w = new double[logWeights.length];
        for(int i = 0; i < w.length; i++){
            w[i] = Math.exp(logWeights[i]-lse);
        }
        return w;
    }

    public static double effectiveSampleSize(double[] normalizedWeights){
        double sumSq = 0.0;
        for (double w : normalizedWeights){
            sumSq += w*w;
        }
        return 1.0/sumSq;
    }

    // matching nodes and weights for Simpson's rule
    // \int f \approx \sum_k weights[k]*f(points[k])
    // numPoints is shifted up by one if it is even, since we
    // need an even number of intervals (and odd number points)
    public static double[][] simpsonNodesAndWeights(int numPoints, double lower, double upper) {
        if (numPoints % 2 == 0) {
            numPoints += 1;
        }
        int n = numPoints - 1; // number of intervals, even
        double h = (upper-lower)/n;

        double[] points = new double[numPoints];
        double[] weights = new double[numPoints];
        for (int i = 0; i <= n; i++) {
            points[i] = lower + i*h;
        }
        weights[0] = h/3.0;
        weights[n] = h/3.0;
        for (int i = 1; i < n; i++) {
            weights[i] = (i%2 == 1) ? (4.0*h/3.0):(2.0*h/3.0);
        }
        return new double[][]{points, weights};
    }

}