package dr.evomodel.treedatalikelihood.continuous;

import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.continuous.ou.canonical.CanonicalNativeBranchGradientCapability;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalBranchMessageContribution;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalLocalTransitionAdjoints;
import dr.evomodel.treedatalikelihood.preorder.BranchSufficientStatistics;
import dr.inference.model.AbstractBlockDiagonalTwoByTwoMatrixParameter;
import dr.inference.model.Parameter;

import java.util.Arrays;

final class OUCanonicalNativeSelectionGradient {

    private final OUCanonicalBranchWiring branchWiring;
    private final OUProcessModel processModel;
    private final CanonicalDebugOptions debugOptions;
    private final CanonicalGradientFallbackPolicy fallbackPolicy;
    private final OUCanonicalBranchFiniteDifference finiteDifference;
    private final OUCanonicalBranchGradientScratch scratch;
    private final int dimension;
    private final ThreadLocal<Integer> debugNodeContext;

    OUCanonicalNativeSelectionGradient(final OUCanonicalBranchWiring branchWiring,
                                       final OUProcessModel processModel,
                                       final CanonicalDebugOptions debugOptions,
                                       final CanonicalGradientFallbackPolicy fallbackPolicy,
                                       final OUCanonicalBranchFiniteDifference finiteDifference,
                                       final OUCanonicalBranchGradientScratch scratch,
                                       final int dimension,
                                       final ThreadLocal<Integer> debugNodeContext) {
        this.branchWiring = branchWiring;
        this.processModel = processModel;
        this.debugOptions = debugOptions;
        this.fallbackPolicy = fallbackPolicy;
        this.finiteDifference = finiteDifference;
        this.scratch = scratch;
        this.dimension = dimension;
        this.debugNodeContext = debugNodeContext;
    }

    double[] getNativeSpecializedGradientWrtSelection(
            final double branchLength,
            final double[] optimum,
            final BranchSufficientStatistics statistics,
            final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter,
            final Parameter requestedParameter) {
        if (!(processModel.getSelectionMatrixParameterization()
                instanceof CanonicalNativeBranchGradientCapability)) {
            throw new IllegalStateException(
                    "Native specialized gradient requires a selection parametrization with native branch capability");
        }

        System.arraycopy(optimum, 0, scratch.stationaryMean, 0, dimension);
        final CanonicalLocalTransitionAdjoints localAdjoints = scratch.localAdjoints;
        branchWiring.fillLocalAdjoints(branchLength, scratch.stationaryMean, statistics, localAdjoints);
        final CanonicalBranchMessageContribution localContribution =
                branchWiring.prepareCanonicalContribution(branchLength, scratch.stationaryMean, statistics);
        final CanonicalNativeBranchGradientCapability nativeBranchGradient =
                (CanonicalNativeBranchGradientCapability) processModel.getSelectionMatrixParameterization();

        OUCanonicalBranchGradientUtils.zero(scratch.compressedNativeGradient);
        OUCanonicalBranchGradientUtils.zero(scratch.rotationGradient);
        final CanonicalLocalTransitionAdjoints adjointsForNative;
        if (fallbackPolicy.transposeNativeDF()) {
            final CanonicalLocalTransitionAdjoints transposed = new CanonicalLocalTransitionAdjoints(dimension);
            OUCanonicalBranchGradientUtils.copyAdjoints(localAdjoints, transposed, true, true);
            OUCanonicalBranchGradientUtils.transposeFlatSquare(localAdjoints.dLogL_dF, transposed.dLogL_dF, dimension);
            adjointsForNative = transposed;
        } else {
            adjointsForNative = localAdjoints;
        }
        try {
            if (fallbackPolicy.useNativeGradientFromContribution()) {
                nativeBranchGradient.accumulateNativeGradientFromCanonicalContribution(
                        processModel.getDiffusionMatrix(),
                        scratch.stationaryMean,
                        branchLength,
                        localContribution,
                        adjointsForNative,
                        scratch.compressedNativeGradient,
                        scratch.rotationGradient);
            } else {
                nativeBranchGradient.accumulateNativeGradientFromAdjoints(
                        processModel.getDiffusionMatrix(),
                        scratch.stationaryMean,
                        branchLength,
                        adjointsForNative,
                        scratch.compressedNativeGradient,
                        scratch.rotationGradient);
            }
        } catch (final RuntimeException e) {
            branchWiring.fillLocalAdjoints(statistics, localAdjoints);
            return finiteDifference.numericalLocalSelectionGradient(
                    branchLength,
                    scratch.stationaryMean,
                    localAdjoints,
                    requestedParameter,
                    e);
        }

        if (!OUCanonicalBranchGradientUtils.isFinite(scratch.compressedNativeGradient)) {
            branchWiring.fillLocalAdjoints(statistics, localAdjoints);
            return finiteDifference.numericalLocalSelectionGradient(
                    branchLength,
                    scratch.stationaryMean,
                    localAdjoints,
                    requestedParameter,
                    new IllegalStateException("Non-finite compressed native gradient"));
        }

        final int nativeDim = blockParameter.getBlockDiagonalNParameters();
        final double[] nativeGradient = new double[nativeDim];
        blockParameter.chainGradient(scratch.compressedNativeGradient, nativeGradient);
        maybeEmitNativeAssemblyDebug(
                statistics,
                branchLength,
                blockParameter,
                requestedParameter,
                scratch.compressedNativeGradient,
                nativeGradient,
                scratch.rotationGradient);
        final double[] direct = CanonicalSelectionGradientProjector.assembleBlockGradientResult(
                dimension, requestedParameter, blockParameter, nativeGradient, scratch.rotationGradient);
        if (debugOptions.isBranchLocalSelectionFiniteDifferenceEnabled()) {
            emitNativeLocalSelectionDebug(branchLength, requestedParameter, direct, localAdjoints);
        }
        if (OUCanonicalBranchGradientUtils.isFinite(direct)) {
            return direct;
        }
        branchWiring.fillLocalAdjoints(statistics, localAdjoints);
        return finiteDifference.numericalLocalSelectionGradient(
                branchLength,
                scratch.stationaryMean,
                localAdjoints,
                requestedParameter,
                new IllegalStateException(
                        "Native specialized tree selection gradient became non-finite: compressed="
                                + Arrays.toString(scratch.compressedNativeGradient)
                                + ", direct=" + Arrays.toString(direct)));
    }

    private void emitNativeLocalSelectionDebug(final double branchLength,
                                               final Parameter requestedParameter,
                                               final double[] direct,
                                               final CanonicalLocalTransitionAdjoints localAdjoints) {
        final double[] numeric = finiteDifference.numericalLocalSelectionGradientGeneric(
                branchLength,
                scratch.stationaryMean,
                localAdjoints,
                requestedParameter);
        double maxAbs = 0.0;
        int maxIdx = -1;
        for (int i = 0; i < Math.min(direct.length, numeric.length); ++i) {
            final double diff = Math.abs(direct[i] - numeric[i]);
            if (diff > maxAbs) {
                maxAbs = diff;
                maxIdx = i;
            }
        }
        System.err.println("NATIVE_SPECIALIZED_LOCAL_FD branchLength=" + branchLength
                + " requested=" + requestedParameter.getId()
                + " maxAbsDiff=" + maxAbs
                + " maxIdx=" + maxIdx
                + " analytic=" + Arrays.toString(direct)
                + " numeric=" + Arrays.toString(numeric));
    }

    private void maybeEmitNativeAssemblyDebug(final BranchSufficientStatistics statistics,
                                              final double branchLength,
                                              final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter,
                                              final Parameter requestedParameter,
                                              final double[] compressedNative,
                                              final double[] nativeBlock,
                                              final double[][] rotation) {
        if (!debugOptions.isNativeAssemblyEnabled()) {
            return;
        }
        final Integer requestedNode = debugOptions.getNativeAssemblyNode();
        if (requestedNode != null) {
            final Integer contextNode = debugNodeContext.get();
            if (contextNode == null || !contextNode.equals(requestedNode)) {
                return;
            }
        }
        final double[] assembled = CanonicalSelectionGradientProjector.assembleBlockGradientResult(
                dimension, requestedParameter, blockParameter, nativeBlock, rotation);
        final Integer contextNode = debugNodeContext.get();
        System.err.println("branchNativeAssemblyDebug node=" + (contextNode == null ? "<unknown>" : contextNode)
                + " branch="
                + statistics.getBranch()
                + " branchLength=" + branchLength
                + " requested=" + requestedParameter.getId()
                + " compressedNative=" + Arrays.toString(compressedNative)
                + " nativeBlock=" + Arrays.toString(nativeBlock)
                + " assembledRequested=" + Arrays.toString(assembled));
    }
}
