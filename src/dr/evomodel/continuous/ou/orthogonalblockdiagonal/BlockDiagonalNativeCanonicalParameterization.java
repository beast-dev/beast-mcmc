package dr.evomodel.continuous.ou.orthogonalblockdiagonal;

import dr.evomodel.continuous.ou.canonical.SpecializedCanonicalSelectionParameterization;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalBranchMessageContribution;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalLocalTransitionAdjoints;
import dr.inference.model.MatrixParameterInterface;

/**
 * Native canonical OU backend for block-diagonal selection parametrizations.
 */
public interface BlockDiagonalNativeCanonicalParameterization
        extends BlockDiagonalCanonicalParameterization, SpecializedCanonicalSelectionParameterization {

    void accumulateNativeGradientFromCanonicalContributionFlat(
            MatrixParameterInterface diffusionMatrix,
            double[] stationaryMean,
            double dt,
            CanonicalBranchMessageContribution contribution,
            CanonicalLocalTransitionAdjoints localAdjoints,
            double[] compressedDAccumulator,
            double[] rotationAccumulator);

    void accumulateNativeGradientFromAdjointsFlat(MatrixParameterInterface diffusionMatrix,
                                                  double[] stationaryMean,
                                                  double dt,
                                                  CanonicalLocalTransitionAdjoints localAdjoints,
                                                  double[] compressedDAccumulator,
                                                  double[] rotationAccumulator);

    void accumulateNativeGradientFromTransitionFlat(double dt,
                                                    double[] stationaryMean,
                                                    double[] dLogL_dF,
                                                    double[] dLogL_df,
                                                    double[] compressedDAccumulator,
                                                    double[] rotationAccumulator);

    void accumulateNativeGradientFromCovarianceStationaryFlat(MatrixParameterInterface diffusionMatrix,
                                                              double dt,
                                                              double[] dLogL_dV,
                                                              double[] compressedDAccumulator,
                                                              double[] rotationAccumulator);
}
