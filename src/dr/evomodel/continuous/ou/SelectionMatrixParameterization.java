package dr.evomodel.continuous.ou;

import dr.inference.model.MatrixParameterInterface;

/**
 * Parametrization-specific access to the OU selection matrix A and its transition map.
 */
public interface SelectionMatrixParameterization {

    MatrixParameterInterface getMatrixParameter();

    int getDimension();

    void fillSelectionMatrixFlat(double[] out);

    void fillTransitionMatrixFlat(double dt, double[] out);

    void accumulateGradientFromTransitionFlat(double dt,
                                              double[] stationaryMean,
                                              double[] dLogL_dF,
                                              double[] dLogL_df,
                                              double[] gradientAccumulator);
}
