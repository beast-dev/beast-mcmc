package dr.evomodel.continuous.ou;

import dr.inference.model.MatrixParameterInterface;

/**
 * Parametrization-specific access to the OU selection matrix A and its transition map.
 */
public interface SelectionMatrixParameterization {

    MatrixParameterInterface getMatrixParameter();

    int getDimension();

    void fillSelectionMatrix(double[][] out);

    void fillTransitionMatrix(double dt, double[][] out);

    void accumulateGradientFromTransition(double dt,
                                          double[] stationaryMean,
                                          double[][] dLogL_dF,
                                          double[] dLogL_df,
                                          double[] gradientAccumulator);
}
