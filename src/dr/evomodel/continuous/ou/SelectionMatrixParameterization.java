package dr.evomodel.continuous.ou;

import dr.inference.model.MatrixParameterInterface;

/**
 * Parametrization-specific access to the OU selection matrix A and its transition map.
 */
public interface SelectionMatrixParameterization {

    MatrixParameterInterface getMatrixParameter();

    int getDimension();

    void fillSelectionMatrix(double[][] out);

    default void fillSelectionMatrixFlat(final double[] out) {
        final int dimension = getDimension();
        final double[][] dense = new double[dimension][dimension];
        fillSelectionMatrix(dense);
        for (int i = 0; i < dimension; ++i) {
            System.arraycopy(dense[i], 0, out, i * dimension, dimension);
        }
    }

    void fillTransitionMatrix(double dt, double[][] out);

    void fillTransitionMatrixFlat(double dt, double[] out);

    void accumulateGradientFromTransition(double dt,
                                          double[] stationaryMean,
                                          double[][] dLogL_dF,
                                          double[] dLogL_df,
                                          double[] gradientAccumulator);

    void accumulateGradientFromTransitionFlat(double dt,
                                              double[] stationaryMean,
                                              double[] dLogL_dF,
                                              double[] dLogL_df,
                                              double[] gradientAccumulator);
}
