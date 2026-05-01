package dr.evomodel.continuous.ou.canonical;

import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalBranchMessageContribution;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalLocalTransitionAdjoints;
import dr.inference.model.MatrixParameterInterface;

/**
 * Capability for branch-local native gradient pullbacks outside the prepared
 * branch traversal path.
 */
public interface CanonicalNativeBranchGradientCapability {

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

    void accumulateNativeGradientFromCanonicalContributionFlat(
            MatrixParameterInterface diffusionMatrix,
            double[] stationaryMean,
            double dt,
            CanonicalBranchMessageContribution contribution,
            CanonicalLocalTransitionAdjoints localAdjoints,
            double[] compressedDAccumulator,
            double[] rotationAccumulator);

    void accumulateNativeGradientFromAdjoints(MatrixParameterInterface diffusionMatrix,
                                             double[] stationaryMean,
                                             double dt,
                                             CanonicalLocalTransitionAdjoints localAdjoints,
                                             double[] compressedDAccumulator,
                                             double[][] rotationAccumulator);

    void accumulateNativeGradientFromAdjointsFlat(MatrixParameterInterface diffusionMatrix,
                                                  double[] stationaryMean,
                                                  double dt,
                                                  CanonicalLocalTransitionAdjoints localAdjoints,
                                                  double[] compressedDAccumulator,
                                                  double[] rotationAccumulator);
}
