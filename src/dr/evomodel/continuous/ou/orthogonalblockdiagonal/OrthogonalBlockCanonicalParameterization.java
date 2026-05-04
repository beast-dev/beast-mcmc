package dr.evomodel.continuous.ou.orthogonalblockdiagonal;

import dr.evomodel.continuous.ou.canonical.SpecializedCanonicalSelectionParameterization;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianTransition;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalBranchMessageContribution;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalLocalTransitionAdjoints;
import dr.inference.model.MatrixParameterInterface;

/**
 * Canonical OU capability exposed by orthogonal block-diagonal selection parametrizations.
 */
public interface OrthogonalBlockCanonicalParameterization extends SpecializedCanonicalSelectionParameterization {

    void fillCanonicalTransition(MatrixParameterInterface diffusionMatrix,
                                 double[] stationaryMean,
                                 double dt,
                                 CanonicalGaussianTransition out);

    OrthogonalBlockPreparedBranchBasis createPreparedBranchBasis();

    void prepareBranchBasis(double dt,
                            double[] stationaryMean,
                            OrthogonalBlockPreparedBranchBasis prepared);

    OrthogonalBlockBranchGradientWorkspace createBranchGradientWorkspace();

    void fillCanonicalTransitionPrepared(OrthogonalBlockPreparedBranchBasis prepared,
                                         MatrixParameterInterface diffusionMatrix,
                                         OrthogonalBlockBranchGradientWorkspace workspace,
                                         CanonicalGaussianTransition out);

    void accumulateNativeGradientFromAdjointsPreparedFlat(
            OrthogonalBlockPreparedBranchBasis prepared,
            MatrixParameterInterface diffusionMatrix,
            CanonicalLocalTransitionAdjoints localAdjoints,
            OrthogonalBlockBranchGradientWorkspace workspace,
            double[] compressedDAccumulator,
            double[] rotationAccumulator);

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

    void accumulateMeanGradient(double dt,
                                double[] dLogL_df,
                                double[] gradientAccumulator);

    double accumulateScalarMeanGradient(double dt,
                                        double[] dLogL_df);

    void accumulateMeanGradientPrepared(OrthogonalBlockPreparedBranchBasis prepared,
                                        double[] dLogL_df,
                                        double[] gradientAccumulator,
                                        OrthogonalBlockBranchGradientWorkspace workspace);

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

    void accumulateDiffusionGradientPreparedFlat(OrthogonalBlockPreparedBranchBasis prepared,
                                                 double[] dLogL_dV,
                                                 boolean transposeAdjoint,
                                                 double[] gradientAccumulator,
                                                 OrthogonalBlockBranchGradientWorkspace workspace);
}
