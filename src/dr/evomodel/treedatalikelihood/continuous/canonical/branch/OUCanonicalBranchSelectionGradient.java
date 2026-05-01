package dr.evomodel.treedatalikelihood.continuous.canonical.branch;

import dr.evomodel.treedatalikelihood.continuous.CanonicalDebugOptions;
import dr.evomodel.treedatalikelihood.continuous.CanonicalGradientFallbackPolicy;
import dr.evomodel.treedatalikelihood.continuous.OUGaussianBranchTransitionProvider;

import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.continuous.ou.canonical.CanonicalNativeBranchGradientCapability;
import dr.evomodel.treedatalikelihood.continuous.canonical.gradient.CanonicalSelectionPullback;
import dr.evomodel.treedatalikelihood.continuous.canonical.gradient.DenseBlockCanonicalSelectionPullback;
import dr.evomodel.treedatalikelihood.continuous.canonical.gradient.DenseCanonicalSelectionPullback;
import dr.evomodel.treedatalikelihood.continuous.canonical.gradient.NativeSpecializedCanonicalSelectionPullback;
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
    private final CanonicalSelectionPullback.DenseSelectionGradientProvider denseGradientProvider;

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
        this.denseGradientProvider = new CanonicalSelectionPullback.DenseSelectionGradientProvider() {
            @Override
            public double[] gradient(final double branchLength,
                                     final double[] optimum,
                                     final BranchSufficientStatistics statistics) {
                return getDenseGradientWrtSelection(branchLength, optimum, statistics);
            }
        };
    }

    void accumulateGlobalRemainderSelectionGradientForBranch(
            final double branchLength,
            final double[] barVdiFlatIn,
            final AbstractBlockDiagonalTwoByTwoMatrixParameter nativeBlockParameter,
            final Parameter requestedParameter,
            final double[] gradientAccumulator) {

        refreshProcessSnapshots();

        final double[] denseGradient = new double[dimension * dimension];
        processModel.accumulateSelectionGradientFromCovarianceFlat(branchLength, barVdiFlatIn, false, denseGradient);
        addInto(selectionPullback(nativeBlockParameter, requestedParameter, false).projectDenseGradient(denseGradient),
                gradientAccumulator);
    }

    double[] getGradientWrtSelection(final double branchLength,
                                     final double[] optimum,
                                     final BranchSufficientStatistics statistics) {
        if (useNumericLocalSelectionGradient()) {
            return finiteDifference.numericalLocalSelectionGradientFromFrozenFactor(
                    branchLength, optimum, statistics, scratch.localAdjoints);
        }
        computeDenseSelectionGradientComponents(branchLength, optimum, statistics);
        if (debugOptions.isSelectionComponentsEnabled()) {
            dr.evomodel.treedatalikelihood.continuous.canonical.CanonicalDiagnosticsLog.warning("selectionComponentsDebug branchLength=" + branchLength
                    + " transition=" + Arrays.toString(scratch.transitionSelectionGradient)
                    + " covariance=" + Arrays.toString(scratch.covarianceSelectionGradient)
                    + " total=" + Arrays.toString(scratch.totalSelectionGradient));
        }
        return scratch.totalSelectionGradient.clone();
    }

    double[] getGradientWrtSelectionComponents(final double branchLength,
                                               final double[] optimum,
                                               final BranchSufficientStatistics statistics) {
        computeDenseSelectionGradientComponents(branchLength, optimum, statistics);
        fillComponentGradient();
        return scratch.selectionComponentGradient.clone();
    }

    double[] getNumericalDenseGradientWrtSelectionComponents(final double branchLength,
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
        final double[] components = new double[3 * totalNumeric.length];
        System.arraycopy(transitionNumeric, 0, components, 0, transitionNumeric.length);
        System.arraycopy(covarianceNumeric, 0, components, transitionNumeric.length, covarianceNumeric.length);
        System.arraycopy(totalNumeric, 0, components, 2 * transitionNumeric.length, totalNumeric.length);
        return components;
    }

    double[] getGradientWrtSelection(final double branchLength,
                                     final double[] optimum,
                                     final BranchSufficientStatistics statistics,
                                     final boolean externalBranch,
                                     final AbstractBlockDiagonalTwoByTwoMatrixParameter nativeBlockParameter,
                                     final Parameter requestedParameter) {
        final boolean useNumeric = useNumericLocalSelectionGradient();
        final CanonicalSelectionPullback pullback = selectionPullback(
                nativeBlockParameter,
                requestedParameter,
                !useNumeric);
        if (useNumeric) {
            return pullback.projectDenseGradient(
                    finiteDifference.numericalLocalSelectionGradientFromFrozenFactor(
                            branchLength, optimum, statistics, scratch.localAdjoints));
        }
        if (debugOptions.isNonOrthogonalOmegaEnabled()
                && nativeBlockParameter != null
                && !(processModel.getSelectionMatrixParameterization()
                instanceof CanonicalNativeBranchGradientCapability)) {
            emitNonNativeBlockDebug(
                branchLength,
                optimum,
                statistics,
                nativeBlockParameter,
                requestedParameter);
        }
        return pullback.gradientForBranch(branchLength, optimum, statistics, denseGradientProvider);
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
        return new DenseCanonicalSelectionPullback(dimension, requestedParameter).projectDenseGradient(
                finiteDifference.numericalLocalSelectionGradientFromFrozenFactor(
                        branchLength, optimum, statistics, scratch.localAdjoints));
    }

    private void computeDenseSelectionGradientComponents(final double branchLength,
                                                         final double[] optimum,
                                                         final BranchSufficientStatistics statistics) {
        refreshProcessSnapshots();
        final CanonicalLocalTransitionAdjoints localAdjoints = scratch.localAdjoints;
        branchWiring.fillLocalAdjoints(branchLength, optimum, statistics, localAdjoints);
        final double[] transitionGradient = scratch.transitionSelectionGradient;
        final double[] covarianceGradient = scratch.covarianceSelectionGradient;
        final double[] totalGradient = scratch.totalSelectionGradient;
        OUCanonicalBranchGradientUtils.zero(transitionGradient);
        OUCanonicalBranchGradientUtils.zero(covarianceGradient);
        OUCanonicalBranchGradientUtils.zero(totalGradient);

        processModel.accumulateSelectionGradientFlat(
                branchLength,
                optimum,
                localAdjoints.dLogL_dF,
                localAdjoints.dLogL_df,
                transitionGradient);
        processModel.accumulateSelectionGradientFromCovarianceFlat(
                branchLength,
                localAdjoints.dLogL_dOmega,
                !fallbackPolicy.useNoTransposeDOmega(),
                covarianceGradient);
        for (int i = 0; i < totalGradient.length; ++i) {
            totalGradient[i] = transitionGradient[i] + covarianceGradient[i];
        }
        if (debugOptions.isBranchLocalSelectionFiniteDifferenceEnabled()) {
            finiteDifference.emitDenseLocalSelectionDebug(branchLength, optimum, localAdjoints, totalGradient);
        }
    }

    private double[] getDenseGradientWrtSelection(final double branchLength,
                                                  final double[] optimum,
                                                  final BranchSufficientStatistics statistics) {
        if (useNumericLocalSelectionGradient()) {
            return finiteDifference.numericalLocalSelectionGradientFromFrozenFactor(
                    branchLength, optimum, statistics, scratch.localAdjoints);
        }
        return getGradientWrtSelection(branchLength, optimum, statistics);
    }

    private CanonicalSelectionPullback selectionPullback(
            final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter,
            final Parameter requestedParameter,
            final boolean allowNativeSpecialized) {
        if (blockParameter == null) {
            return new DenseCanonicalSelectionPullback(dimension, requestedParameter);
        }
        if (allowNativeSpecialized
                && processModel.getSelectionMatrixParameterization()
                instanceof CanonicalNativeBranchGradientCapability
                && !fallbackPolicy.useDensePullbackForOrthogonal()) {
            return new NativeSpecializedCanonicalSelectionPullback(
                    nativeSelectionGradient,
                    blockParameter,
                    requestedParameter);
        }
        return new DenseBlockCanonicalSelectionPullback(
                dimension,
                requestedParameter,
                blockParameter);
    }

    private void emitNonNativeBlockDebug(final double branchLength,
                                         final double[] optimum,
                                         final BranchSufficientStatistics statistics,
                                         final AbstractBlockDiagonalTwoByTwoMatrixParameter nativeBlockParameter,
                                         final Parameter requestedParameter) {
        if (debugOptions.isNonOrthogonalOmegaEnabled()) {
            branchWiring.fillLocalAdjoints(statistics, scratch.localAdjoints);
            dr.evomodel.treedatalikelihood.continuous.canonical.CanonicalDiagnosticsLog.warning("NONORTH OMEGA branchLength=" + branchLength
                    + " dOmega=" + Arrays.toString(scratch.localAdjoints.dLogL_dOmega)
                    + " drift=" + Arrays.toString(processModel.getDriftMatrix().getParameterAsMatrix()[0]));
        }
        final double[] denseGradient = getGradientWrtSelection(branchLength, optimum, statistics);
        final double[] analytic = selectionPullback(
                nativeBlockParameter,
                requestedParameter,
                false).projectDenseGradient(denseGradient);
        if (debugOptions.isNonOrthogonalOmegaEnabled()) {
            branchWiring.fillLocalAdjoints(statistics, scratch.localAdjoints);
            final double[] numeric = finiteDifference.numericalLocalSelectionGradientGeneric(
                    branchLength,
                    optimum,
                    scratch.localAdjoints,
                    requestedParameter);
            dr.evomodel.treedatalikelihood.continuous.canonical.CanonicalDiagnosticsLog.warning("NONORTH LOCAL branchLength=" + branchLength
                    + " analytic=" + Arrays.toString(analytic)
                    + " numeric=" + Arrays.toString(numeric));
        }
    }

    private boolean useNumericLocalSelectionGradient() {
        return fallbackPolicy.useNumericLocalSelectionGradient();
    }

    private void refreshProcessSnapshots() {
        branchTransitionProvider.getProcessModel();
    }

    private void fillComponentGradient() {
        final int blockLength = scratch.transitionSelectionGradient.length;
        System.arraycopy(scratch.transitionSelectionGradient, 0, scratch.selectionComponentGradient, 0, blockLength);
        System.arraycopy(scratch.covarianceSelectionGradient, 0, scratch.selectionComponentGradient, blockLength, blockLength);
        System.arraycopy(scratch.totalSelectionGradient, 0, scratch.selectionComponentGradient, 2 * blockLength, blockLength);
    }

    private static void addInto(final double[] values, final double[] accumulator) {
        for (int k = 0; k < accumulator.length; ++k) {
            accumulator[k] += values[k];
        }
    }
}
