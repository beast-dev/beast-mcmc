package dr.evomodel.continuous.ou.orthogonalblockdiagonal;

import dr.evomodel.continuous.ou.SpecializedCanonicalSelectionParameterization;
import dr.evomodel.treedatalikelihood.continuous.gaussian.CanonicalGaussianTransition;
import dr.evomodel.treedatalikelihood.continuous.gaussian.message.CanonicalBranchMessageContribution;
import dr.evomodel.treedatalikelihood.continuous.gaussian.message.CanonicalLocalTransitionAdjoints;
import dr.inference.model.MatrixParameterInterface;

/**
 * Canonical OU capability exposed by orthogonal block-diagonal selection parametrizations.
 */
public interface OrthogonalBlockCanonicalParameterization extends SpecializedCanonicalSelectionParameterization {

    void fillTransitionCovariance(MatrixParameterInterface diffusionMatrix,
                                  double dt,
                                  double[][] out);

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

    void accumulateNativeGradientFromTransition(double dt,
                                                double[] stationaryMean,
                                                double[][] dLogL_dF,
                                                double[] dLogL_df,
                                                double[] compressedDAccumulator,
                                                double[][] rotationAccumulator);

    void accumulateNativeGradientFromCovarianceStationary(MatrixParameterInterface diffusionMatrix,
                                                          double dt,
                                                          double[][] dLogL_dV,
                                                          double[] compressedDAccumulator,
                                                          double[][] rotationAccumulator);

    void accumulateDiffusionGradient(MatrixParameterInterface diffusionMatrix,
                                     double dt,
                                     double[][] dLogL_dV,
                                     double[] gradientAccumulator);

    void accumulateNativeGradientFromCanonicalContribution(
            MatrixParameterInterface diffusionMatrix,
            double[] stationaryMean,
            double dt,
            CanonicalBranchMessageContribution contribution,
            CanonicalLocalTransitionAdjoints localAdjoints,
            double[] compressedDAccumulator,
            double[][] rotationAccumulator);

    void accumulateNativeGradientFromAdjoints(MatrixParameterInterface diffusionMatrix,
                                             double[] stationaryMean,
                                             double dt,
                                             CanonicalLocalTransitionAdjoints localAdjoints,
                                             double[] compressedDAccumulator,
                                             double[][] rotationAccumulator);

    void accumulateMeanGradient(double dt,
                                double[] dLogL_df,
                                double[] gradientAccumulator);

    double accumulateScalarMeanGradient(double dt,
                                        double[] dLogL_df);

    void accumulateMeanGradientPrepared(OrthogonalBlockPreparedBranchBasis prepared,
                                        double[] dLogL_df,
                                        double[] gradientAccumulator,
                                        OrthogonalBlockBranchGradientWorkspace workspace);

    void accumulateDiffusionGradientPreparedFlat(OrthogonalBlockPreparedBranchBasis prepared,
                                                 double[] dLogL_dV,
                                                 boolean transposeAdjoint,
                                                 double[] gradientAccumulator,
                                                 OrthogonalBlockBranchGradientWorkspace workspace);
}
