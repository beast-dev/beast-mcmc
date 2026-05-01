package dr.evomodel.treedatalikelihood.continuous;

import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.continuous.ou.canonical.CanonicalNativeBranchGradientCapability;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalLocalTransitionAdjoints;
import dr.evomodel.treedatalikelihood.preorder.BranchSufficientStatistics;
import dr.inference.model.AbstractBlockDiagonalTwoByTwoMatrixParameter;
import dr.inference.model.Parameter;

import java.util.Arrays;

final class OUCanonicalBranchSelectionGradient {

    private final OUGaussianBranchTransitionProvider branchTransitionProvider;
    private final OUCanonicalBranchWiring branchWiring;
    private final OUProcessModel processModel;
    private final CanonicalDebugOptions debugOptions;
    private final CanonicalGradientFallbackPolicy fallbackPolicy;
    private final OUCanonicalBranchFiniteDifference finiteDifference;
    private final OUCanonicalNativeSelectionGradient nativeSelectionGradient;
    private final OUCanonicalBranchGradientScratch scratch;
    private final int dimension;

    OUCanonicalBranchSelectionGradient(
            final OUGaussianBranchTransitionProvider branchTransitionProvider,
            final OUCanonicalBranchWiring branchWiring,
            final OUProcessModel processModel,
            final CanonicalDebugOptions debugOptions,
            final CanonicalGradientFallbackPolicy fallbackPolicy,
            final OUCanonicalBranchFiniteDifference finiteDifference,
            final OUCanonicalNativeSelectionGradient nativeSelectionGradient,
            final OUCanonicalBranchGradientScratch scratch,
            final int dimension) {
        this.branchTransitionProvider = branchTransitionProvider;
        this.branchWiring = branchWiring;
        this.processModel = processModel;
        this.debugOptions = debugOptions;
        this.fallbackPolicy = fallbackPolicy;
        this.finiteDifference = finiteDifference;
        this.nativeSelectionGradient = nativeSelectionGradient;
        this.scratch = scratch;
        this.dimension = dimension;
    }

    void accumulateGlobalRemainderSelectionGradientForBranch(
            final double branchLength,
            final double[] barVdiFlatIn,
            final AbstractBlockDiagonalTwoByTwoMatrixParameter nativeBlockParameter,
            final Parameter requestedParameter,
            final double[] gradientAccumulator) {

        refreshProcessSnapshots();

        final double[][] barVdi2D = new double[dimension][dimension];
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                barVdi2D[i][j] = 0.5 * (
                        barVdiFlatIn[i * dimension + j] + barVdiFlatIn[j * dimension + i]);
            }
        }

        if (nativeBlockParameter == null) {
            if (requestedParameter == null) {
                processModel.accumulateSelectionGradientFromCovariance(
                        branchLength, barVdi2D, gradientAccumulator);
            } else {
                final double[] denseGradient = new double[dimension * dimension];
                processModel.accumulateSelectionGradientFromCovariance(branchLength, barVdi2D, denseGradient);
                final double[] reordered = CanonicalSelectionGradientProjector.reorderDenseGradientForRequestedParameter(
                        dimension, requestedParameter, denseGradient);
                addInto(reordered, gradientAccumulator);
            }
            return;
        }

        final double[] denseGradient = new double[dimension * dimension];
        processModel.accumulateSelectionGradientFromCovariance(branchLength, barVdi2D, denseGradient);
        final double[] pulled = CanonicalSelectionGradientProjector.pullBackDenseGradientToBlock(
                dimension, requestedParameter, nativeBlockParameter, denseGradient);
        addInto(pulled, gradientAccumulator);
    }

    double[] getGradientWrtSelection(final double branchLength,
                                     final double[] optimum,
                                     final BranchSufficientStatistics statistics) {
        if (useNumericLocalSelectionGradient()) {
            return finiteDifference.numericalLocalSelectionGradientFromFrozenFactor(
                    branchLength, optimum, statistics, scratch.localAdjoints);
        }
        final double[][] components = computeDenseSelectionGradientComponents(branchLength, optimum, statistics);
        if (debugOptions.isSelectionComponentsEnabled()) {
            System.err.println("selectionComponentsDebug branchLength=" + branchLength
                    + " transition=" + Arrays.toString(components[0])
                    + " covariance=" + Arrays.toString(components[1])
                    + " total=" + Arrays.toString(components[2]));
        }
        return components[2].clone();
    }

    double[][] getGradientWrtSelectionComponents(final double branchLength,
                                                 final double[] optimum,
                                                 final BranchSufficientStatistics statistics) {
        final double[][] components = computeDenseSelectionGradientComponents(branchLength, optimum, statistics);
        return new double[][]{
                components[0].clone(),
                components[1].clone(),
                components[2].clone()
        };
    }

    double[][] getNumericalDenseGradientWrtSelectionComponents(final double branchLength,
                                                               final double[] optimum,
                                                               final BranchSufficientStatistics statistics) {
        refreshProcessSnapshots();
        final CanonicalLocalTransitionAdjoints localAdjoints = scratch.localAdjoints;
        branchWiring.fillLocalAdjoints(branchLength, optimum, statistics, localAdjoints);
        final CanonicalLocalTransitionAdjoints transitionOnly = new CanonicalLocalTransitionAdjoints(dimension);
        final CanonicalLocalTransitionAdjoints covarianceOnly = new CanonicalLocalTransitionAdjoints(dimension);
        OUCanonicalBranchGradientUtils.copyAdjoints(localAdjoints, transitionOnly, true, false);
        OUCanonicalBranchGradientUtils.copyAdjoints(localAdjoints, covarianceOnly, false, true);
        final double[] transitionNumeric = finiteDifference.numericalLocalSelectionGradientDenseA(
                branchLength, optimum, transitionOnly);
        final double[] covarianceNumeric = finiteDifference.numericalLocalSelectionGradientDenseA(
                branchLength, optimum, covarianceOnly);
        final double[] totalNumeric = new double[transitionNumeric.length];
        for (int i = 0; i < totalNumeric.length; ++i) {
            totalNumeric[i] = transitionNumeric[i] + covarianceNumeric[i];
        }
        return new double[][]{transitionNumeric, covarianceNumeric, totalNumeric};
    }

    double[] getGradientWrtSelection(final double branchLength,
                                     final double[] optimum,
                                     final BranchSufficientStatistics statistics,
                                     final boolean externalBranch,
                                     final AbstractBlockDiagonalTwoByTwoMatrixParameter nativeBlockParameter,
                                     final Parameter requestedParameter) {
        if (useNumericLocalSelectionGradient()) {
            final double[] denseGradient = finiteDifference.numericalLocalSelectionGradientFromFrozenFactor(
                    branchLength, optimum, statistics, scratch.localAdjoints);
            if (nativeBlockParameter == null) {
                return CanonicalSelectionGradientProjector.reorderDenseGradientForRequestedParameter(
                        dimension, requestedParameter, denseGradient);
            }
            return CanonicalSelectionGradientProjector.pullBackDenseGradientToBlock(
                    dimension, requestedParameter, nativeBlockParameter, denseGradient);
        }
        if (nativeBlockParameter == null) {
            final double[] denseGradient = getGradientWrtSelection(branchLength, optimum, statistics);
            return CanonicalSelectionGradientProjector.reorderDenseGradientForRequestedParameter(
                    dimension, requestedParameter, denseGradient);
        }

        if (!(processModel.getSelectionMatrixParameterization()
                instanceof CanonicalNativeBranchGradientCapability)) {
            return getNonNativeBlockGradient(branchLength, optimum, statistics, nativeBlockParameter, requestedParameter);
        }

        if (fallbackPolicy.useDensePullbackForOrthogonal()) {
            return getGradientWrtSelectionViaDensePullback(
                    branchLength, optimum, statistics, nativeBlockParameter, requestedParameter);
        }

        return getNativeSpecializedGradientWrtSelection(
                branchLength,
                optimum,
                statistics,
                nativeBlockParameter,
                requestedParameter);
    }

    double[] getNumericalFrozenLocalGradientWrtSelection(
            final double branchLength,
            final double[] optimum,
            final BranchSufficientStatistics statistics,
            final AbstractBlockDiagonalTwoByTwoMatrixParameter nativeBlockParameter,
            final Parameter requestedParameter) {
        if (nativeBlockParameter != null) {
            if (fallbackPolicy.useNumericLocalSelectionFromFrozenLogFactor()) {
                return finiteDifference.numericalLocalSelectionGradientGenericFromFrozenFactor(
                        branchLength, optimum, statistics, requestedParameter);
            }
            branchWiring.fillLocalAdjoints(branchLength, optimum, statistics, scratch.localAdjoints);
            return finiteDifference.numericalLocalSelectionGradientGeneric(
                    branchLength,
                    optimum,
                    scratch.localAdjoints,
                    requestedParameter);
        }
        return CanonicalSelectionGradientProjector.reorderDenseGradientForRequestedParameter(
                dimension,
                requestedParameter,
                finiteDifference.numericalLocalSelectionGradientFromFrozenFactor(
                        branchLength, optimum, statistics, scratch.localAdjoints));
    }

    double[] getGradientWrtSelectionViaDensePullback(
            final double branchLength,
            final double[] optimum,
            final BranchSufficientStatistics statistics,
            final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter,
            final Parameter requestedParameter) {
        final double[] denseGradient = getGradientWrtSelection(branchLength, optimum, statistics);
        if (blockParameter == null) {
            return CanonicalSelectionGradientProjector.reorderDenseGradientForRequestedParameter(
                    dimension, requestedParameter, denseGradient);
        }
        return CanonicalSelectionGradientProjector.pullBackDenseGradientToBlock(
                dimension, requestedParameter, blockParameter, denseGradient);
    }

    double[] getNativeSpecializedGradientWrtSelection(
            final double branchLength,
            final double[] optimum,
            final BranchSufficientStatistics statistics,
            final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter,
            final Parameter requestedParameter) {
        return nativeSelectionGradient.getNativeSpecializedGradientWrtSelection(
                branchLength,
                optimum,
                statistics,
                blockParameter,
                requestedParameter);
    }

    private double[][] computeDenseSelectionGradientComponents(final double branchLength,
                                                               final double[] optimum,
                                                               final BranchSufficientStatistics statistics) {
        refreshProcessSnapshots();
        final CanonicalLocalTransitionAdjoints localAdjoints = scratch.localAdjoints;
        branchWiring.fillLocalAdjoints(branchLength, optimum, statistics, localAdjoints);
        final double[] transitionGradient = new double[scratch.matrixGradient.length];
        final double[] covarianceGradient = new double[scratch.matrixGradient.length];
        final double[] totalGradient = new double[scratch.matrixGradient.length];

        OUCanonicalBranchGradientUtils.copyFromFlatInto(
                localAdjoints.dLogL_dF, scratch.transitionAdjoint, dimension);
        processModel.accumulateSelectionGradient(
                branchLength,
                optimum,
                scratch.transitionAdjoint,
                localAdjoints.dLogL_df,
                transitionGradient);
        if (fallbackPolicy.useNoTransposeDOmega()) {
            OUCanonicalBranchGradientUtils.copyFromFlatInto(
                    localAdjoints.dLogL_dOmega, scratch.covarianceAdjoint, dimension);
        } else {
            OUCanonicalBranchGradientUtils.transposeFromFlatInto(
                    localAdjoints.dLogL_dOmega, scratch.covarianceAdjoint, dimension);
        }
        processModel.accumulateSelectionGradientFromCovariance(
                branchLength,
                scratch.covarianceAdjoint,
                covarianceGradient);
        for (int i = 0; i < totalGradient.length; ++i) {
            totalGradient[i] = transitionGradient[i] + covarianceGradient[i];
        }
        if (debugOptions.isBranchLocalSelectionFiniteDifferenceEnabled()) {
            finiteDifference.emitDenseLocalSelectionDebug(branchLength, optimum, localAdjoints, totalGradient);
        }
        return new double[][]{transitionGradient, covarianceGradient, totalGradient};
    }

    private double[] getNonNativeBlockGradient(final double branchLength,
                                               final double[] optimum,
                                               final BranchSufficientStatistics statistics,
                                               final AbstractBlockDiagonalTwoByTwoMatrixParameter nativeBlockParameter,
                                               final Parameter requestedParameter) {
        if (debugOptions.isNonOrthogonalOmegaEnabled()) {
            branchWiring.fillLocalAdjoints(statistics, scratch.localAdjoints);
            System.err.println("NONORTH OMEGA branchLength=" + branchLength
                    + " dOmega=" + Arrays.toString(scratch.localAdjoints.dLogL_dOmega)
                    + " drift=" + Arrays.toString(processModel.getDriftMatrix().getParameterAsMatrix()[0]));
        }
        final double[] denseGradient = getGradientWrtSelection(branchLength, optimum, statistics);
        final double[] analytic = CanonicalSelectionGradientProjector.pullBackDenseGradientToBlock(
                dimension, requestedParameter, nativeBlockParameter, denseGradient);
        if (debugOptions.isNonOrthogonalOmegaEnabled()) {
            branchWiring.fillLocalAdjoints(statistics, scratch.localAdjoints);
            final double[] numeric = finiteDifference.numericalLocalSelectionGradientGeneric(
                    branchLength,
                    optimum,
                    scratch.localAdjoints,
                    requestedParameter);
            System.err.println("NONORTH LOCAL branchLength=" + branchLength
                    + " analytic=" + Arrays.toString(analytic)
                    + " numeric=" + Arrays.toString(numeric));
        }
        return analytic;
    }

    private boolean useNumericLocalSelectionGradient() {
        return fallbackPolicy.useNumericLocalSelectionGradient();
    }

    private void refreshProcessSnapshots() {
        branchTransitionProvider.getProcessModel();
    }

    private static void addInto(final double[] values, final double[] accumulator) {
        for (int k = 0; k < accumulator.length; ++k) {
            accumulator[k] += values[k];
        }
    }
}
