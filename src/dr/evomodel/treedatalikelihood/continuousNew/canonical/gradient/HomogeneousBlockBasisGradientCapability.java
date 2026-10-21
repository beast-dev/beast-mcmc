package dr.evomodel.treedatalikelihood.continuousNew.canonical.gradient;

import dr.evomodel.continuous.ou.blockdiagonal.BlockBasisCanonicalAdjointBuilder;
import dr.evomodel.continuous.ou.blockdiagonal.BlockBasisGlobalSnapshot;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.message.CanonicalLocalTransitionAdjoints;
import dr.inference.model.AbstractBlockDiagonalTwoByTwoMatrixParameter;

/**
 * Capability exposed by homogeneous block-basis OU transition providers.
 */
public interface HomogeneousBlockBasisGradientCapability {

    BlockBasisGlobalSnapshot getBlockBasisSnapshot();

    BlockBasisCanonicalAdjointBuilder getBlockBasisAdjointBuilder();

    AbstractBlockDiagonalTwoByTwoMatrixParameter getBlockParameter();

    double contractBranchLengthGradient(double branchLength,
                                        CanonicalLocalTransitionAdjoints localAdjoints);
}
