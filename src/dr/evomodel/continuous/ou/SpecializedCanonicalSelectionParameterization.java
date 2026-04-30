package dr.evomodel.continuous.ou;

import dr.evomodel.treedatalikelihood.continuous.gaussian.CanonicalGaussianTransition;
import dr.evomodel.treedatalikelihood.continuous.gaussian.message.CanonicalLocalTransitionAdjoints;
import dr.inference.model.MatrixParameterInterface;

/**
 * Optional canonical OU backend for selection parametrizations with a faster
 * native transition, covariance, or gradient path.
 */
public interface SpecializedCanonicalSelectionParameterization extends SelectionMatrixParameterization {

    int getSelectionGradientDimension();

    int getCompressedSelectionGradientDimension();

    int getNativeSelectionGradientScratchDimension();

    void finishNativeSelectionGradient(double[] compressedGradient,
                                       double[] nativeGradientScratch,
                                       double[] rotationGradientFlat,
                                       double[] gradientOut);

    void fillTransitionCovariance(MatrixParameterInterface diffusionMatrix,
                                  double dt,
                                  double[][] out);

    void fillCanonicalTransition(MatrixParameterInterface diffusionMatrix,
                                 double[] stationaryMean,
                                 double dt,
                                 CanonicalGaussianTransition out);

    CanonicalPreparedBranchHandle createPreparedBranchHandle();

    CanonicalBranchWorkspace createBranchWorkspace();

    void prepareBranch(double dt,
                       double[] stationaryMean,
                       CanonicalPreparedBranchHandle prepared);

    void fillCanonicalTransitionPrepared(CanonicalPreparedBranchHandle prepared,
                                         MatrixParameterInterface diffusionMatrix,
                                         CanonicalBranchWorkspace workspace,
                                         CanonicalGaussianTransition out);

    void accumulateNativeGradientFromAdjointsPreparedFlat(
            CanonicalPreparedBranchHandle prepared,
            MatrixParameterInterface diffusionMatrix,
            CanonicalLocalTransitionAdjoints localAdjoints,
            CanonicalBranchWorkspace workspace,
            double[] compressedDAccumulator,
            double[] rotationAccumulator);

    void accumulateMeanGradientPrepared(CanonicalPreparedBranchHandle prepared,
                                        double[] dLogL_df,
                                        double[] gradientAccumulator,
                                        CanonicalBranchWorkspace workspace);

    void accumulateDiffusionGradientPreparedFlat(CanonicalPreparedBranchHandle prepared,
                                                 double[] dLogL_dV,
                                                 boolean transposeAdjoint,
                                                 double[] gradientAccumulator,
                                                 CanonicalBranchWorkspace workspace);
}
