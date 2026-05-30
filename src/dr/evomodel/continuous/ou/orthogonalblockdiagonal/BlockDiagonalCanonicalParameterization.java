package dr.evomodel.continuous.ou.orthogonalblockdiagonal;

import dr.evomodel.continuous.ou.SelectionMatrixParameterization;
import dr.evomodel.continuous.ou.canonical.CanonicalPreparedTransitionCapability;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalBranchMessageContribution;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianTransition;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalLocalTransitionAdjoints;
import dr.inference.model.MatrixParameterInterface;

/**
 * Canonical OU capability exposed by block-diagonal selection parametrizations.
 */
public interface BlockDiagonalCanonicalParameterization
        extends SelectionMatrixParameterization, CanonicalPreparedTransitionCapability {

    void fillCanonicalTransition(MatrixParameterInterface diffusionMatrix,
                                 double[] stationaryMean,
                                 double dt,
                                 CanonicalGaussianTransition out);

    BlockDiagonalPreparedBranchBasis createPreparedBranchBasis();

    void prepareBranchBasis(double dt,
                            double[] stationaryMean,
                            BlockDiagonalPreparedBranchBasis prepared);

    BlockDiagonalBranchGradientWorkspace createBranchGradientWorkspace();

    void fillCanonicalTransitionPrepared(BlockDiagonalPreparedBranchBasis prepared,
                                         MatrixParameterInterface diffusionMatrix,
                                         BlockDiagonalBranchGradientWorkspace workspace,
                                         CanonicalGaussianTransition out);

    void fillCanonicalLocalAdjoints(MatrixParameterInterface diffusionMatrix,
                                    double[] stationaryMean,
                                    double dt,
                                    CanonicalBranchMessageContribution contribution,
                                    CanonicalLocalTransitionAdjoints out);

    void accumulateMeanGradient(double dt,
                                double[] dLogL_df,
                                double[] gradientAccumulator);

    double accumulateScalarMeanGradient(double dt,
                                        double[] dLogL_df);
}
