package dr.math.distributions;

import dr.math.MathUtils;
import java.util.Arrays;

public class PointMassMixtureDistribution implements MultivariateDistribution, RandomGenerator {

    public static final String TYPE = "PointMassMixture";

    //
    private final double[] weights;
    // Each column of realizedValues matrix is a vector representing multivariate parameter value
    private final double[][] realizedValues;
    boolean weightsNormalized;
    private double[] probabilites;


    public PointMassMixtureDistribution(double[] weights, double[][] realizedValues, boolean weightsNormalized) {
        this.weights = weights;
        this.realizedValues = realizedValues;
        this.weightsNormalized = weightsNormalized;
        this.probabilites = computeProbabilities(weights,weightsNormalized);
    }

    public String getType() {
        return TYPE;
    }

    public double[] computeProbabilities(double[] weights, boolean weightsNormalized){
        if(weightsNormalized){
            return weights;
        }else{
            double normConst = 0;
            for(int i = 0; i < weights.length; i++){
                normConst = normConst + weights[i];
            }
            for(int j = 0; j < weights.length; j++){
                weights[j] = weights[j]/normConst;
            }
            return weights;
        }
    }

    public double logPdf(double[] x){
        return logPdf(x,probabilites,realizedValues);
    }

    public static double logPdf(double[] x, double[] probs, double[][] values){
        int numRows = values.length;
        int numCols = values[0].length;

        for(int j = 0; j < numCols; j++){
            double[] parameterValue = new double[numRows];
            for(int i = 0; i < numRows; i++){
                parameterValue[i] = values[i][j];
                //System.err.println("parameter " + j + " coordinate " + i + " is: " + parameterValue[i]);
                //System.err.println("x coord " + x[i]);
            }
            //System.err.println("probability parameter " + j + " is: " + probs[j]);
            //System.err.println("x: " + x[0] + " paramVal: " + parameterValue[0] + " x.equals(ParameterValue): " + Arrays.equals(x,parameterValue));
            if(Arrays.equals(x,parameterValue)) {
                return Math.log(probs[j]);
            }
        }
        // x is not equal to any of realized values
        return Math.log(0);
    }

    public double[][] getCovarianceMatrix(){
        return getCovarianceMatrix(probabilites,realizedValues);
    }

    public static double[][] getCovarianceMatrix(double[] probs, double[][] values){

        int numRows = values.length;
        int numCols = values[0].length;
        double[][] cov = new double[numRows][numCols];

        for(int i = 0; i < numRows; i++){
            for(int j = 0; j < numCols; j++){
                cov = addMatrices(cov,
                        multMatrixByScalar(probs[i]*probs[j],
                                multColumnVecByRowVec(getCol(values,i),getCol(values,j))));
            }
        }
        cov = addMatrices(cov,
                multMatrixByScalar(-1,
                        multColumnVecByRowVec(getMean(probs, values),getMean(probs, values))));

        return cov;
    }

    public double[][] getScaleMatrix(){
        return getScaleMatrix(probabilites,realizedValues);
    }

    public static double[][] getScaleMatrix(double[] probs, double[][] values){
        return getCovarianceMatrix(probs, values);
    }

    public double[] getMean(){
        return getMean(probabilites,realizedValues);
    }

    public static double[] getMean(double[] probs, double[][] values){
        double[] mean = new double[values.length];
        for(int i = 0; i < values.length; i++){
            for(int j = 0; j < values[i].length; j++){
                mean[i] = mean[i] + probs[j]*values[i][j];
            }
        }
        return mean;
    }



    public static double[][] multMatrixByScalar(double scalar, double[][] mat){
        for(int i = 0; i < mat.length; i++){
            for(int j =0; j < mat[i].length; j++){
                mat[i][j] = scalar*mat[i][j];
            }
        }
        return mat;
    }

    public static double[][] multColumnVecByRowVec(double[] rowVec, double[] columnVec){
        double[][] prod = new double[rowVec.length][columnVec.length];
        for(int i = 0; i < rowVec.length; i++){
            for(int j = 0; j < columnVec.length; j++){
                prod[i][j] = columnVec[i]*rowVec[j];
            }
        }
        return prod;
    }

    public static double[][] addMatrices(double[][] matA, double[][] matB){
        int numRows = matA.length;
        int numCols = matA[0].length;
        for(int i = 0; i < numRows; i++){
            for(int j = 0; j < numCols; j++){
                matA[i][j] = matA[i][j] + matB[i][j];
            }
        }
        return matA;
    }

    public static double[] getCol(double[][] mat, int index){
        double[] col = new double[mat.length];
        for(int i = 0; i < mat.length; i++){
            col[i] = mat[i][index];
        }
        return col;
    }


    // RandomGenerator interface
    public Object nextRandom() {
        //for(int i = 0; i < 4; i++) {
        //    System.err.println("probabilities[" + i + "]: " + probabilites[i]);
        //}

        int index = MathUtils.randomChoicePDF(probabilites);

        //System.err.println("index: " + index);

        //for(int i = 0; i < 6; i++) {
        //    System.err.println("realizedValues[" + i + "]: " + realizedValues[i][index]);
        //}

        return getCol(realizedValues, index);
    }

    public double logPdf(Object x) {
        double[] v = (double[]) x;
        return logPdf(v);
    }

}