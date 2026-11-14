package dr.evomodel.treedatalikelihood.continuousNew.canonical.gradient;

import dr.evomodel.treedatalikelihood.continuousNew.canonical.message.CanonicalLocalTransitionAdjoints;

/**
 * Capability exposed by homogeneous block-basis OU transition providers.
 *
 * <p>The gradient layer stays independent of the block-diagonal backend: the
 * provider, which owns the backend snapshot and adjoint builder, creates the
 * selection-gradient pullback itself.
 */
public interface HomogeneousBlockBasisGradientCapability {

    CanonicalSelectionGradientPullback createBlockBasisSelectionGradientPullback(int dimension);

    double contractBranchLengthGradient(double branchLength,
                                        CanonicalLocalTransitionAdjoints localAdjoints);
}
