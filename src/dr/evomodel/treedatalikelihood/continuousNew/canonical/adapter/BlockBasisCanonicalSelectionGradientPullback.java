package dr.evomodel.treedatalikelihood.continuousNew.canonical.adapter;

import dr.evomodel.continuous.canonical.math.MatrixOps;
import dr.evomodel.continuous.canonical.workspace.BranchAdjointWorkspace;
import dr.evomodel.continuous.canonical.workspace.BranchGradientWorkspace;
import dr.evomodel.continuous.canonical.workspace.GradientPullbackWorkspace;
import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.continuous.ou.blockdiagonal.BlockBasisCanonicalAdjointBuilder;
import dr.evomodel.continuous.ou.blockdiagonal.BlockBasisGlobalSnapshot;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.gradient.BranchGradientInputs;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.gradient.CanonicalSelectionGradientProjector;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.gradient.CanonicalSelectionGradientPullback;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.message.CanonicalGaussianMessageOps;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.message.CanonicalLocalTransitionAdjoints;
import dr.inference.model.AbstractBlockDiagonalTwoByTwoMatrixParameter;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import java.util.Arrays;

/**
 * Pulls homogeneous block-basis OU adjoints back to native block, basis, and
 * diffusion parameters.
 */
public final class BlockBasisCanonicalSelectionGradientPullback implements CanonicalSelectionGradientPullback {

    private final HomogeneousBlockBasisCanonicalOUBranchTransitionProvider blockBasis;
    private final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter;
    private final BlockBasisCanonicalAdjointBuilder adjointBuilder;
    private final int dimension;
    private final double[] compressedDGradient;
    private final double[] nativeBlockGradient;
    private final DenseMatrix64F stationaryCovarianceGradient;
    private final DenseMatrix64F rotationGradient;
    private final double[] rootMeanGradientDBasis;
    private int observedTipCount;

    BlockBasisCanonicalSelectionGradientPullback(
            final HomogeneousBlockBasisCanonicalOUBranchTransitionProvider blockBasis,
            final int dimension) {
        this.blockBasis = blockBasis;
        this.blockParameter = blockBasis.getBlockParameter();
        this.adjointBuilder = blockBasis.getBlockBasisAdjointBuilder();
        this.dimension = dimension;
        if (blockParameter.getRowDimension() != dimension
                || blockParameter.getColumnDimension() != dimension) {
            throw new IllegalArgumentException("Block-basis gradient dimension mismatch.");
        }
        this.compressedDGradient = new double[blockParameter.getCompressedDDimension()];
        this.nativeBlockGradient = new double[blockParameter.getBlockDiagonalNParameters()];
        this.stationaryCovarianceGradient = new DenseMatrix64F(dimension, dimension);
        this.rotationGradient = new DenseMatrix64F(dimension, dimension);
        this.rootMeanGradientDBasis = new double[dimension];
    }

    @Override
    public void initialize(final BranchGradientWorkspace workspace,
                           final double[] gradA,
                           final double[] gradMu) {
        blockBasis.getBlockBasisSnapshot();
        Arrays.fill(compressedDGradient, 0.0);
        Arrays.fill(nativeBlockGradient, 0.0);
        CommonOps.fill(stationaryCovarianceGradient, 0.0);
        CommonOps.fill(rotationGradient, 0.0);
        Arrays.fill(rootMeanGradientDBasis, 0.0);
        observedTipCount = 0;
        if (gradMu.length != dimension) {
            throw new IllegalArgumentException(
                    "Block-basis stationary-mean gradient length must be " + dimension
                            + ", found " + gradMu.length);
        }
    }

    @Override
    public void clearWorkerBuffers(final BranchGradientWorkspace workspace,
                                   final int gradALength,
                                   final int gradMuLength) {
        throw new UnsupportedOperationException(
                "Homogeneous block-basis OU gradients currently run through the sequential reducer.");
    }

    @Override
    public void prepareWorkspace(final BranchGradientWorkspace workspace) { }

    @Override
    public void accumulateForBranch(final OUProcessModel processModel,
                                    final BranchGradientInputs inputs,
                                    final int activeIndex,
                                    final CanonicalLocalTransitionAdjoints localAdjoints,
                                    final BranchGradientWorkspace workspace,
                                    final double[] gradA,
                                    final double[] gradQ,
                                    final double[] gradMu) {
        final BlockBasisGlobalSnapshot snapshot = blockBasis.getBlockBasisSnapshot();
        adjointBuilder.accumulateBranchPullback(
                snapshot,
                inputs.getBranchLength(activeIndex),
                localAdjoints,
                compressedDGradient,
                stationaryCovarianceGradient);

        if (localAdjoints.hasObservedTipValueAdjoint) {
            adjointBuilder.addTipTransformRGradient(
                    snapshot,
                    localAdjoints.observedTipValue,
                    localAdjoints.dLogL_dObservedTipValue,
                    rotationGradient);
            accumulateStationaryMeanTransformGradient(
                    snapshot,
                    localAdjoints.dLogL_dObservedTipValue,
                    -1.0,
                    gradMu);
            observedTipCount++;
        }

        if (localAdjoints.hasObservationBasisAdjoint) {
            addFlatToDense(localAdjoints.dLogL_dObservationBasis, rotationGradient);
            for (int i = 0; i < dimension; ++i) {
                gradMu[i] += localAdjoints.dLogL_dStationaryMeanDirect[i];
            }
        }
    }

    @Override
    public void reduceWorker(final BranchGradientWorkspace worker,
                             final BranchGradientWorkspace reductionWorkspace,
                             final double[] gradA) {
        throw new UnsupportedOperationException(
                "Homogeneous block-basis OU gradients currently run through the sequential reducer.");
    }

    @Override
    public void finish(final BranchGradientInputs inputs,
                       final BranchGradientWorkspace workspace,
                       final double[] gradA,
                       final double[] gradQ,
                       final double[] gradMu) {
        final BlockBasisGlobalSnapshot snapshot = blockBasis.getBlockBasisSnapshot();
        if (inputs.getRootDiffusionScale() > 0.0) {
            accumulateRootPriorPullback(inputs, workspace, gradQ, gradMu);
        }

        adjointBuilder.finalizeGlobalPullback(
                snapshot,
                stationaryCovarianceGradient,
                compressedDGradient,
                rotationGradient,
                gradQ);
        adjointBuilder.addJacobianTermGradient(snapshot, observedTipCount, rotationGradient);
        blockParameter.chainGradient(compressedDGradient, nativeBlockGradient);

        final double[] assembledGradient =
                CanonicalSelectionGradientProjector.assembleBlockGradientResultFlat(
                        dimension,
                        blockParameter.getParameter(),
                        blockParameter,
                        nativeBlockGradient,
                        rotationGradient.data);
        if (assembledGradient.length != gradA.length) {
            throw new IllegalArgumentException(
                    "Block-basis selection gradient length " + assembledGradient.length
                            + " does not match requested length " + gradA.length + ".");
        }
        System.arraycopy(assembledGradient, 0, gradA, 0, assembledGradient.length);
    }

    @Override
    public boolean handlesRootDiffusionGradient() {
        return true;
    }

    private void accumulateRootPriorPullback(final BranchGradientInputs inputs,
                                             final BranchGradientWorkspace workspace,
                                             final double[] gradQ,
                                             final double[] gradMu) {
        final BranchAdjointWorkspace adjoint = workspace.adjoint;
        final GradientPullbackWorkspace gradient = workspace.gradient;

        CanonicalGaussianMessageOps.combineStates(
                inputs.getRootPreOrderState(),
                inputs.getRootPostOrderState(),
                adjoint.combinedState);

        adjoint.fillMomentsFromCanonical(
                adjoint.combinedState,
                adjoint.mean,
                adjoint.covarianceFlat,
                dimension);
        adjoint.fillMomentsFromCanonical(
                inputs.getRootPreOrderState(),
                adjoint.mean2,
                gradient.covariance2,
                dimension);

        final double[] priorPrecision = inputs.getRootPreOrderState().precision;
        for (int i = 0; i < dimension; ++i) {
            double rootGradient = 0.0;
            final int iOff = i * dimension;
            for (int j = 0; j < dimension; ++j) {
                rootGradient += priorPrecision[iOff + j] * (adjoint.mean[j] - adjoint.mean2[j]);
            }
            rootMeanGradientDBasis[i] = rootGradient;
        }
        final BlockBasisGlobalSnapshot snapshot = blockBasis.getBlockBasisSnapshot();
        adjointBuilder.addTipTransformRGradient(
                snapshot,
                adjoint.mean2,
                rootMeanGradientDBasis,
                rotationGradient);
        accumulateStationaryMeanTransformGradient(
                snapshot,
                rootMeanGradientDBasis,
                -1.0,
                gradMu);

        for (int i = 0; i < dimension; ++i) {
            final double deltaI = adjoint.mean[i] - adjoint.mean2[i];
            final int iOff = i * dimension;
            for (int j = 0; j < dimension; ++j) {
                gradient.covarianceAdjointFlat[iOff + j] =
                        adjoint.covarianceFlat[iOff + j]
                                + deltaI * (adjoint.mean[j] - adjoint.mean2[j]);
            }
        }

        MatrixOps.matMul(
                priorPrecision,
                gradient.covarianceAdjointFlat,
                gradient.matrixProductFlat,
                dimension);
        MatrixOps.matMul(
                gradient.matrixProductFlat,
                priorPrecision,
                gradient.covarianceAdjointFlat,
                dimension);

        final double rootScale = inputs.getRootDiffusionScale();
        for (int i = 0; i < dimension; ++i) {
            final int iOff = i * dimension;
            for (int j = 0; j < dimension; ++j) {
                final int ij = iOff + j;
                gradient.covarianceAdjointFlat[ij] =
                        rootScale * (-0.5 * priorPrecision[ij]
                                + 0.5 * gradient.covarianceAdjointFlat[ij]);
            }
        }
        adjointBuilder.accumulateTransformedDiffusionPullback(
                snapshot,
                gradient.covarianceAdjointFlat,
                rotationGradient,
                gradQ);
    }

    private void accumulateStationaryMeanTransformGradient(final BlockBasisGlobalSnapshot snapshot,
                                                           final double[] gradientDBasis,
                                                           final double scale,
                                                           final double[] gradMu) {
        final DenseMatrix64F t = snapshot.getRotationInverse();
        for (int j = 0; j < dimension; ++j) {
            double sum = 0.0;
            for (int i = 0; i < dimension; ++i) {
                sum += t.get(i, j) * gradientDBasis[i];
            }
            gradMu[j] += scale * sum;
        }
    }

    private void addFlatToDense(final double[] source, final DenseMatrix64F target) {
        if (source.length < dimension * dimension) {
            throw new IllegalArgumentException("Flat matrix gradient is too short.");
        }
        for (int i = 0; i < dimension * dimension; ++i) {
            target.data[i] += source[i];
        }
    }
}
