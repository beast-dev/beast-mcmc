package dr.math.distributions;

import java.util.Arrays;

public class WishartSufficientStatistics {

    public WishartSufficientStatistics(int dim) {
        df = 0;
        scaleMatrix = new double[dim][dim];
    }

    public WishartSufficientStatistics(int df, double[][] scaleMatrix) {
        this.df = df;
        this.scaleMatrix = scaleMatrix;
    }

    public final int getDf() {
        return df;
    }

    public final double[][] getScaleMatrix() {
        return scaleMatrix;
    }
    
    public final void incrementDf(int n) {
        df += n;
    }

    public final void clear() {
        df = 0;
        for (double[] v : scaleMatrix) {
            Arrays.fill(v, 0.0);
        }
    }

    private int df;
    private final double[][] scaleMatrix;
}