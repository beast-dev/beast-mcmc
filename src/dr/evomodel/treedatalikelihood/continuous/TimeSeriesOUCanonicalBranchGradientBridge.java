package dr.evomodel.treedatalikelihood.continuous;

import dr.evomodel.treedatalikelihood.preorder.BranchSufficientStatistics;
import dr.evomodel.treedatalikelihood.preorder.MatrixSufficientStatistics;
import dr.inference.model.AbstractBlockDiagonalTwoByTwoMatrixParameter;
import dr.inference.model.OrthogonalMatrixProvider;
import dr.inference.model.Parameter;
import dr.inference.model.TransposedMatrixParameter;
import dr.evomodel.treedatalikelihood.continuous.gaussian.message.CanonicalBranchMessageContribution;
import dr.evomodel.treedatalikelihood.continuous.gaussian.message.CanonicalLocalTransitionAdjoints;
import dr.evomodel.continuous.ou.orthogonalblockdiagonal.OrthogonalBlockCanonicalParameterization;
import dr.evomodel.continuous.ou.OUProcessModel;
import org.ejml.data.DenseMatrix64F;

import java.util.Arrays;

/**
 * Tree-side bridge that reconstructs a branch-local canonical posterior from
 * the preorder/postorder sufficient statistics and then reuses the time-series
 * canonical adjoint machinery for OU branch gradients.
 */
public final class TimeSeriesOUCanonicalBranchGradientBridge {
    private static final ThreadLocal<Integer> DEBUG_NODE_CONTEXT = new ThreadLocal<Integer>();

    private final TimeSeriesOUGaussianBranchTransitionProvider branchTransitionProvider;
    private final OUProcessModel processModel;
    private final int dimension;
    private final OUCanonicalBranchWiring branchWiring;
    private final CanonicalDebugOptions debugOptions;
    private final CanonicalGradientFallbackPolicy fallbackPolicy;

    private final CanonicalLocalTransitionAdjoints localAdjoints;

    private final double[] branchMeanGradient;
    private final double[] matrixGradient;
    private final double[] stationaryMean;
    private final double[] compressedNativeGradient;
    private final double[][] rotationGradient;
    private final double[][] transitionAdjointScratch;
    private final double[][] covarianceAdjointScratch;
    private final double[] parameterRestoreScratch;

    public TimeSeriesOUCanonicalBranchGradientBridge(
            final TimeSeriesOUGaussianBranchTransitionProvider branchTransitionProvider) {
        this(branchTransitionProvider,
                CanonicalDebugOptions.fromSystemProperties(),
                CanonicalGradientFallbackPolicy.fromSystemProperties());
    }

    TimeSeriesOUCanonicalBranchGradientBridge(
            final TimeSeriesOUGaussianBranchTransitionProvider branchTransitionProvider,
            final CanonicalDebugOptions debugOptions,
            final CanonicalGradientFallbackPolicy fallbackPolicy) {
        if (branchTransitionProvider == null) {
            throw new IllegalArgumentException("branchTransitionProvider must not be null");
        }
        if (debugOptions == null) {
            throw new IllegalArgumentException("debugOptions must not be null");
        }
        if (fallbackPolicy == null) {
            throw new IllegalArgumentException("fallbackPolicy must not be null");
        }
        this.branchTransitionProvider = branchTransitionProvider;
        this.debugOptions = debugOptions;
        this.fallbackPolicy = fallbackPolicy;
        this.branchWiring = new OUCanonicalBranchWiring(branchTransitionProvider, debugOptions, fallbackPolicy);
        this.processModel = branchWiring.getProcessModel();
        this.dimension = branchWiring.getDimension();
        this.localAdjoints = new CanonicalLocalTransitionAdjoints(dimension);
        this.branchMeanGradient = new double[dimension];
        this.matrixGradient = new double[dimension * dimension];
        this.stationaryMean = new double[dimension];
        this.compressedNativeGradient = new double[dimension + 2 * (dimension / 2)];
        this.rotationGradient = new double[dimension][dimension];
        this.transitionAdjointScratch = new double[dimension][dimension];
        this.covarianceAdjointScratch = new double[dimension][dimension];
        this.parameterRestoreScratch = new double[Math.max(1, dimension * dimension)];
    }

    public static void setDebugNodeContext(final Integer nodeNumber) {
        if (nodeNumber == null) {
            DEBUG_NODE_CONTEXT.remove();
        } else {
            DEBUG_NODE_CONTEXT.set(nodeNumber);
        }
    }

    public int getDimension() {
        return dimension;
    }

    public void fillLocalAdjoints(final BranchSufficientStatistics statistics,
                                  final CanonicalLocalTransitionAdjoints out) {
        refreshProcessSnapshots();
        branchWiring.fillLocalAdjoints(statistics, out);
    }

    public double[] getGradientWrtBranchDrift(final BranchSufficientStatistics statistics) {
        fillLocalAdjoints(statistics, localAdjoints);
        branchWiring.accumulateBranchMeanGradient(statistics.getBranch(), localAdjoints.dLogL_df, branchMeanGradient);
        return branchMeanGradient.clone();
    }

    public double[] getGradientWrtBranchDrift(final double branchLength,
                                              final double[] optimum,
                                              final BranchSufficientStatistics statistics,
                                              final int nodeNumber) {
        fillLocalAdjoints(statistics, localAdjoints);
        branchWiring.accumulateBranchMeanGradient(statistics.getBranch(), localAdjoints.dLogL_df, branchMeanGradient);
        maybeEmitLocalMeanFDDebug(branchLength, optimum, statistics, nodeNumber, branchMeanGradient);
        return branchMeanGradient.clone();
    }

    public double[] getGradientWrtVariance(final double branchLength,
                                           final double[] optimum,
                                           final BranchSufficientStatistics statistics) {
        // Keep process-model diffusion snapshot synchronized with current tree parameters.
        // Without this refresh, accumulateDiffusionGradient can use stale Q while adjoints
        // are computed from current branch statistics.
        refreshProcessSnapshots();
        branchWiring.fillLocalAdjoints(branchLength, optimum, statistics, localAdjoints);
        zero(matrixGradient);
        if (processModel.getSelectionMatrixParameterization()
                instanceof OrthogonalBlockCanonicalParameterization) {
            // Branch-local canonical adjoints in the tree bridge use the opposite
            // covariance orientation relative to the orthogonal exact helper.
            // Match the same orientation convention used by the selection-covariance path.
            transposeFromFlatInto(localAdjoints.dLogL_dOmega, covarianceAdjointScratch, dimension);
            ((OrthogonalBlockCanonicalParameterization)
                    processModel.getSelectionMatrixParameterization())
                    .accumulateDiffusionGradient(
                            processModel.getDiffusionMatrix(),
                            branchLength,
                            covarianceAdjointScratch,
                            matrixGradient);
        } else {
            copyFromFlatInto(localAdjoints.dLogL_dOmega, covarianceAdjointScratch, dimension);
            processModel.accumulateDiffusionGradient(branchLength, covarianceAdjointScratch, matrixGradient);
        }
        return matrixGradient.clone();
    }

    /**
     * Accumulates into {@code gradientAccumulator} the branch contribution from the
     * pruning-remainder covariance adjoint {@code barVdiFlatIn}.
     */
    public void accumulateGlobalRemainderSelectionGradientForBranch(
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

        if (nativeBlockParameter != null
                && processModel.getSelectionMatrixParameterization()
                instanceof OrthogonalBlockCanonicalParameterization) {
            final double[] denseGradient = new double[dimension * dimension];
            processModel.accumulateSelectionGradientFromCovariance(branchLength, barVdi2D, denseGradient);
            final double[] pulled = pullBackDenseGradientToBlock(
                    requestedParameter, nativeBlockParameter, denseGradient);
            for (int k = 0; k < gradientAccumulator.length; ++k) {
                gradientAccumulator[k] += pulled[k];
            }
            return;
        }

        if (nativeBlockParameter == null) {
            if (requestedParameter == null) {
                processModel.accumulateSelectionGradientFromCovariance(
                        branchLength, barVdi2D, gradientAccumulator);
            } else {
                final double[] denseGradient = new double[dimension * dimension];
                processModel.accumulateSelectionGradientFromCovariance(branchLength, barVdi2D, denseGradient);
                final double[] reordered = reorderDenseGradientForRequestedParameter(
                        requestedParameter, denseGradient);
                for (int k = 0; k < gradientAccumulator.length; ++k) {
                    gradientAccumulator[k] += reordered[k];
                }
            }
            return;
        }

        final double[] denseGradient = new double[dimension * dimension];
        processModel.accumulateSelectionGradientFromCovariance(branchLength, barVdi2D, denseGradient);
        final double[] pulled = pullBackDenseGradientToBlock(
                requestedParameter, nativeBlockParameter, denseGradient);
        for (int k = 0; k < gradientAccumulator.length; ++k) {
            gradientAccumulator[k] += pulled[k];
        }
    }

    public double[] getGradientWrtSelection(final double branchLength,
                                            final double[] optimum,
                                            final BranchSufficientStatistics statistics) {
        if (useNumericLocalSelectionGradient()) {
            return numericalLocalSelectionGradientFromFrozenFactor(branchLength, optimum, statistics);
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

    public double[][] getGradientWrtSelectionComponents(final double branchLength,
                                                        final double[] optimum,
                                                        final BranchSufficientStatistics statistics) {
        final double[][] components = computeDenseSelectionGradientComponents(branchLength, optimum, statistics);
        return new double[][]{
                components[0].clone(),
                components[1].clone(),
                components[2].clone()
        };
    }

    public double[][] getNumericalDenseGradientWrtSelectionComponents(final double branchLength,
                                                                      final double[] optimum,
                                                                      final BranchSufficientStatistics statistics) {
        refreshProcessSnapshots();
        branchWiring.fillLocalAdjoints(branchLength, optimum, statistics, localAdjoints);
        final CanonicalLocalTransitionAdjoints transitionOnly = new CanonicalLocalTransitionAdjoints(dimension);
        final CanonicalLocalTransitionAdjoints covarianceOnly = new CanonicalLocalTransitionAdjoints(dimension);
        copyAdjoints(localAdjoints, transitionOnly, true, false);
        copyAdjoints(localAdjoints, covarianceOnly, false, true);
        final double[] transitionNumeric = numericalLocalSelectionGradientDenseA(
                branchLength, optimum, transitionOnly);
        final double[] covarianceNumeric = numericalLocalSelectionGradientDenseA(
                branchLength, optimum, covarianceOnly);
        final double[] totalNumeric = new double[transitionNumeric.length];
        for (int i = 0; i < totalNumeric.length; ++i) {
            totalNumeric[i] = transitionNumeric[i] + covarianceNumeric[i];
        }
        return new double[][]{transitionNumeric, covarianceNumeric, totalNumeric};
    }

    private double[][] computeDenseSelectionGradientComponents(final double branchLength,
                                                               final double[] optimum,
                                                               final BranchSufficientStatistics statistics) {
        refreshProcessSnapshots();
        branchWiring.fillLocalAdjoints(branchLength, optimum, statistics, localAdjoints);
        final double[] transitionGradient = new double[matrixGradient.length];
        final double[] covarianceGradient = new double[matrixGradient.length];
        final double[] totalGradient = new double[matrixGradient.length];

        copyFromFlatInto(localAdjoints.dLogL_dF, transitionAdjointScratch, dimension);
        processModel.accumulateSelectionGradient(
                branchLength,
                optimum,
                transitionAdjointScratch,
                localAdjoints.dLogL_df,
                transitionGradient);
        if (fallbackPolicy.useNoTransposeDOmega()) {
            copyFromFlatInto(localAdjoints.dLogL_dOmega, covarianceAdjointScratch, dimension);
        } else {
            transposeFromFlatInto(localAdjoints.dLogL_dOmega, covarianceAdjointScratch, dimension);
        }
        processModel.accumulateSelectionGradientFromCovariance(
                branchLength,
                covarianceAdjointScratch,
                covarianceGradient);
        for (int i = 0; i < totalGradient.length; ++i) {
            totalGradient[i] = transitionGradient[i] + covarianceGradient[i];
        }
        if (debugOptions.isBranchLocalSelectionFiniteDifferenceEnabled()) {
            emitDenseLocalSelectionDebug(branchLength, optimum, localAdjoints, totalGradient);
        }
        return new double[][]{transitionGradient, covarianceGradient, totalGradient};
    }

    public double[] getGradientWrtSelection(final double branchLength,
                                            final double[] optimum,
                                            final BranchSufficientStatistics statistics,
                                            final boolean externalBranch,
                                            final AbstractBlockDiagonalTwoByTwoMatrixParameter nativeBlockParameter,
                                            final Parameter requestedParameter) {
        if (useNumericLocalSelectionGradient()) {
            final double[] denseGradient = numericalLocalSelectionGradientFromFrozenFactor(branchLength, optimum, statistics);
            if (nativeBlockParameter == null) {
                return reorderDenseGradientForRequestedParameter(requestedParameter, denseGradient);
            }
            return pullBackDenseGradientToBlock(requestedParameter, nativeBlockParameter, denseGradient);
        }
        if (nativeBlockParameter == null) {
            final double[] denseGradient = getGradientWrtSelection(branchLength, optimum, statistics);
            return reorderDenseGradientForRequestedParameter(requestedParameter, denseGradient);
        }

        if (!(processModel.getSelectionMatrixParameterization()
                instanceof OrthogonalBlockCanonicalParameterization)) {
            if (debugOptions.isNonOrthogonalOmegaEnabled()) {
                fillLocalAdjoints(statistics, localAdjoints);
                System.err.println("NONORTH OMEGA branchLength=" + branchLength
                        + " dOmega=" + Arrays.toString(localAdjoints.dLogL_dOmega)
                        + " drift=" + Arrays.toString(processModel.getDriftMatrix().getParameterAsMatrix()[0]));
            }
            final double[] denseGradient = getGradientWrtSelection(branchLength, optimum, statistics);
            final double[] analytic = pullBackDenseGradientToBlock(requestedParameter, nativeBlockParameter, denseGradient);
            if (debugOptions.isNonOrthogonalOmegaEnabled()) {
                fillLocalAdjoints(statistics, localAdjoints);
                final double[] numeric = numericalLocalSelectionGradientGeneric(
                        branchLength,
                        optimum,
                        localAdjoints,
                        requestedParameter);
                System.err.println("NONORTH LOCAL branchLength=" + branchLength
                        + " analytic=" + Arrays.toString(analytic)
                        + " numeric=" + Arrays.toString(numeric));
            }
            return analytic;
        }

        if (fallbackPolicy.useDensePullbackForOrthogonal()) {
            final double[] denseGradient = getGradientWrtSelection(branchLength, optimum, statistics);
            return pullBackDenseGradientToBlock(requestedParameter, nativeBlockParameter, denseGradient);
        }

        return getNativeOrthogonalBlockGradientWrtSelection(
                branchLength,
                optimum,
                statistics,
                nativeBlockParameter,
                requestedParameter);
    }

    public double[] getNumericalFrozenLocalGradientWrtSelection(final double branchLength,
                                                                final double[] optimum,
                                                                final BranchSufficientStatistics statistics,
                                                                final AbstractBlockDiagonalTwoByTwoMatrixParameter nativeBlockParameter,
                                                                final Parameter requestedParameter) {
        if (nativeBlockParameter != null) {
            if (fallbackPolicy.useNumericLocalSelectionFromFrozenLogFactor()) {
                return numericalLocalSelectionGradientGenericFromFrozenFactor(
                        branchLength, optimum, statistics, requestedParameter);
            }
            branchWiring.fillLocalAdjoints(branchLength, optimum, statistics, localAdjoints);
            return numericalLocalSelectionGradientGeneric(
                    branchLength,
                    optimum,
                    localAdjoints,
                    requestedParameter);
        }
        return reorderDenseGradientForRequestedParameter(
                requestedParameter,
                numericalLocalSelectionGradientFromFrozenFactor(branchLength, optimum, statistics));
    }

    public double[] getGradientWrtSelectionViaDensePullback(final double branchLength,
                                                            final double[] optimum,
                                                            final BranchSufficientStatistics statistics,
                                                            final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter,
                                                            final Parameter requestedParameter) {
        final double[] denseGradient = getGradientWrtSelection(branchLength, optimum, statistics);
        if (blockParameter == null) {
            return reorderDenseGradientForRequestedParameter(requestedParameter, denseGradient);
        }
        return pullBackDenseGradientToBlock(requestedParameter, blockParameter, denseGradient);
    }

    public double[] getNativeOrthogonalBlockGradientWrtSelection(final double branchLength,
                                                                 final double[] optimum,
                                                                 final BranchSufficientStatistics statistics,
                                                                 final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter,
                                                                 final Parameter requestedParameter) {
        if (!(processModel.getSelectionMatrixParameterization()
                instanceof OrthogonalBlockCanonicalParameterization)) {
            throw new IllegalStateException("Native orthogonal block gradient requires orthogonal block selection parametrization");
        }

        System.arraycopy(optimum, 0, stationaryMean, 0, dimension);
        // Ensure local adjoints are refreshed for this exact branch contribution before pullback.
        branchWiring.fillLocalAdjoints(branchLength, stationaryMean, statistics, localAdjoints);
        final CanonicalBranchMessageContribution localContribution =
                branchWiring.prepareCanonicalContribution(branchLength, stationaryMean, statistics);
        final OrthogonalBlockCanonicalParameterization orthogonalParameterization =
                (OrthogonalBlockCanonicalParameterization) processModel.getSelectionMatrixParameterization();

        zero(compressedNativeGradient);
        zero(rotationGradient);
        final CanonicalLocalTransitionAdjoints adjointsForNative;
        if (fallbackPolicy.transposeNativeDF()) {
            final CanonicalLocalTransitionAdjoints transposed = new CanonicalLocalTransitionAdjoints(dimension);
            copyAdjoints(localAdjoints, transposed, true, true);
            transposeFlatSquare(localAdjoints.dLogL_dF, transposed.dLogL_dF, dimension);
            adjointsForNative = transposed;
        } else {
            adjointsForNative = localAdjoints;
        }
        try {
            if (fallbackPolicy.useNativeGradientFromContribution()) {
                orthogonalParameterization.accumulateNativeGradientFromCanonicalContribution(
                        processModel.getDiffusionMatrix(),
                        stationaryMean,
                        branchLength,
                        localContribution,
                        adjointsForNative,
                        compressedNativeGradient,
                        rotationGradient);
            } else {
                orthogonalParameterization.accumulateNativeGradientFromAdjoints(
                        processModel.getDiffusionMatrix(),
                        stationaryMean,
                        branchLength,
                        adjointsForNative,
                        compressedNativeGradient,
                        rotationGradient);
            }
        } catch (final RuntimeException e) {
            branchWiring.fillLocalAdjoints(statistics, localAdjoints);
            return numericalLocalSelectionGradient(
                    orthogonalParameterization,
                    branchLength,
                    stationaryMean,
                    localAdjoints,
                    requestedParameter,
                    e);
        }

        if (!isFinite(compressedNativeGradient)) {
            branchWiring.fillLocalAdjoints(statistics, localAdjoints);
            return numericalLocalSelectionGradient(
                    orthogonalParameterization,
                    branchLength,
                    stationaryMean,
                    localAdjoints,
                    requestedParameter,
                    new IllegalStateException("Non-finite compressed native gradient"));
        }

        final int nativeDim = blockParameter.getBlockDiagonalNParameters();
        final double[] nativeGradient = new double[nativeDim];
        blockParameter.chainGradient(compressedNativeGradient, nativeGradient);
        maybeEmitNativeAssemblyDebug(
                statistics,
                branchLength,
                blockParameter,
                requestedParameter,
                compressedNativeGradient,
                nativeGradient,
                rotationGradient);
        final double[] direct = assembleBlockGradientResult(requestedParameter, blockParameter, nativeGradient, rotationGradient);
        if (debugOptions.isBranchLocalSelectionFiniteDifferenceEnabled()) {
            final double[] numeric = numericalLocalSelectionGradientGeneric(
                    branchLength,
                    stationaryMean,
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
            System.err.println("NATIVE_ORTH_LOCAL_FD branchLength=" + branchLength
                    + " requested=" + requestedParameter.getId()
                    + " maxAbsDiff=" + maxAbs
                    + " maxIdx=" + maxIdx
                    + " analytic=" + Arrays.toString(direct)
                    + " numeric=" + Arrays.toString(numeric));
        }
        if (isFinite(direct)) {
            return direct;
        }
        branchWiring.fillLocalAdjoints(statistics, localAdjoints);
        return numericalLocalSelectionGradient(
                orthogonalParameterization,
                branchLength,
                stationaryMean,
                localAdjoints,
                requestedParameter,
                new IllegalStateException(
                        "Native orthogonal block tree selection gradient became non-finite: compressed="
                                + Arrays.toString(compressedNativeGradient)
                + ", direct=" + Arrays.toString(direct)));
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
            final Integer contextNode = DEBUG_NODE_CONTEXT.get();
            if (contextNode == null || !contextNode.equals(requestedNode)) {
                return;
            }
        }
        final double[] assembled = assembleBlockGradientResult(requestedParameter, blockParameter, nativeBlock, rotation);
        final Integer contextNode = DEBUG_NODE_CONTEXT.get();
        System.err.println("branchNativeAssemblyDebug node=" + (contextNode == null ? "<unknown>" : contextNode)
                + " branch="
                + statistics.getBranch()
                + " branchLength=" + branchLength
                + " requested=" + requestedParameter.getId()
                + " compressedNative=" + Arrays.toString(compressedNative)
                + " nativeBlock=" + Arrays.toString(nativeBlock)
                + " assembledRequested=" + Arrays.toString(assembled));
    }

    private double[] numericalLocalSelectionGradient(
            final OrthogonalBlockCanonicalParameterization orthogonalParameterization,
            final double branchLength,
            final double[] optimum,
            final CanonicalLocalTransitionAdjoints adjoints,
            final Parameter requestedParameter,
            final RuntimeException cause) {
        final int parameterDimension = requestedParameter.getDimension();
        ensureParameterScratchCapacity(parameterDimension);
        for (int i = 0; i < parameterDimension; ++i) {
            parameterRestoreScratch[i] = requestedParameter.getParameterValue(i);
        }

        final double[] gradient = new double[parameterDimension];
        try {
            for (int i = 0; i < parameterDimension; ++i) {
                final double base = parameterRestoreScratch[i];
                final double step = 1.0e-6 * Math.max(1.0, Math.abs(base));

                requestedParameter.setParameterValue(i, base + step);
                final double plus = evaluateLocalSelectionScore(branchLength, optimum, adjoints);

                requestedParameter.setParameterValue(i, base - step);
                final double minus = evaluateLocalSelectionScore(branchLength, optimum, adjoints);

                if (!Double.isFinite(plus) || !Double.isFinite(minus)) {
                    throw new IllegalStateException(
                            "Non-finite local selection score in numerical fallback at index=" + i
                                    + ", base=" + base
                                    + ", step=" + step
                                    + ", plus=" + plus
                                    + ", minus=" + minus);
                }

                gradient[i] = (plus - minus) / (2.0 * step);
                if (!Double.isFinite(gradient[i])) {
                    throw new IllegalStateException(
                            "Non-finite local numerical gradient at index=" + i
                                    + ", plus=" + plus
                                    + ", minus=" + minus
                                    + ", step=" + step);
                }
                requestedParameter.setParameterValue(i, base);
            }
            return gradient;
        } catch (final RuntimeException e) {
            for (int i = 0; i < parameterDimension; ++i) {
                requestedParameter.setParameterValue(i, parameterRestoreScratch[i]);
            }
            final IllegalStateException wrapped = new IllegalStateException(
                    "Failed numerical local fallback after analytic native block selection gradient failure",
                    e);
            if (cause != null) {
                wrapped.addSuppressed(cause);
            }
            throw wrapped;
        } finally {
            for (int i = 0; i < parameterDimension; ++i) {
                if (requestedParameter.getParameterValue(i) != parameterRestoreScratch[i]) {
                    requestedParameter.setParameterValue(i, parameterRestoreScratch[i]);
                }
            }
        }
    }

    private double evaluateLocalSelectionScore(
            final double branchLength,
            final double[] optimum,
            final CanonicalLocalTransitionAdjoints adjoints) {
        refreshProcessSnapshots();
        return branchWiring.evaluateLocalSelectionScore(branchLength, optimum, adjoints);
    }

    /**
     * Keep the time-series process snapshots synchronized with current tree-side
     * diffusion/selection parameters before canonical bridge gradient evaluation.
     */
    private void refreshProcessSnapshots() {
        branchTransitionProvider.getProcessModel();
    }

    private double[] numericalLocalSelectionGradientGeneric(
            final double branchLength,
            final double[] optimum,
            final CanonicalLocalTransitionAdjoints adjoints,
            final Parameter requestedParameter) {
        final int parameterDimension = requestedParameter.getDimension();
        final double[] restore = new double[parameterDimension];
        for (int i = 0; i < parameterDimension; ++i) {
            restore[i] = requestedParameter.getParameterValue(i);
        }

        final double[] gradient = new double[parameterDimension];
        try {
            for (int i = 0; i < parameterDimension; ++i) {
                final double base = restore[i];
                final double step = 1.0e-6 * Math.max(1.0, Math.abs(base));

                requestedParameter.setParameterValue(i, base + step);
                final double plus = evaluateLocalSelectionScore(branchLength, optimum, adjoints);

                requestedParameter.setParameterValue(i, base - step);
                final double minus = evaluateLocalSelectionScore(branchLength, optimum, adjoints);

                gradient[i] = (plus - minus) / (2.0 * step);
                requestedParameter.setParameterValue(i, base);
            }
            return gradient;
        } finally {
            for (int i = 0; i < parameterDimension; ++i) {
                if (requestedParameter.getParameterValue(i) != restore[i]) {
                    requestedParameter.setParameterValue(i, restore[i]);
                }
            }
        }
    }

    private double[] numericalLocalSelectionGradientGenericFromFrozenFactor(
            final double branchLength,
            final double[] optimum,
            final BranchSufficientStatistics statistics,
            final Parameter requestedParameter) {
        final int parameterDimension = requestedParameter.getDimension();
        final double[] restore = new double[parameterDimension];
        for (int i = 0; i < parameterDimension; ++i) {
            restore[i] = requestedParameter.getParameterValue(i);
        }

        final double[] gradient = new double[parameterDimension];
        try {
            for (int i = 0; i < parameterDimension; ++i) {
                final double base = restore[i];
                final double step = 1.0e-6 * Math.max(1.0, Math.abs(base));

                requestedParameter.setParameterValue(i, base + step);
                refreshProcessSnapshots();
                final double plus = branchWiring.evaluateFrozenLocalLogFactor(
                        branchLength, optimum, statistics.getAbove(), statistics.getBelow());

                requestedParameter.setParameterValue(i, base - step);
                refreshProcessSnapshots();
                final double minus = branchWiring.evaluateFrozenLocalLogFactor(
                        branchLength, optimum, statistics.getAbove(), statistics.getBelow());

                gradient[i] = (plus - minus) / (2.0 * step);
                requestedParameter.setParameterValue(i, base);
            }
            return gradient;
        } finally {
            for (int i = 0; i < parameterDimension; ++i) {
                if (requestedParameter.getParameterValue(i) != restore[i]) {
                    requestedParameter.setParameterValue(i, restore[i]);
                }
            }
            refreshProcessSnapshots();
        }
    }



    private static boolean isFinite(final double[] values) {
        for (double value : values) {
            if (!Double.isFinite(value)) {
                return false;
            }
        }
        return true;
    }

    private static boolean isFinite(final double[][] values) {
        for (double[] row : values) {
            if (!isFinite(row)) {
                return false;
            }
        }
        return true;
    }

    private static boolean hasNaN(final double[] values) {
        for (double value : values) {
            if (Double.isNaN(value)) {
                return true;
            }
        }
        return false;
    }

    private static boolean hasNaN(final double[][] values) {
        for (double[] row : values) {
            if (hasNaN(row)) {
                return true;
            }
        }
        return false;
    }

    private static void assertNoNaNStatistics(final BranchSufficientStatistics statistics) {
        final MatrixSufficientStatistics branch = statistics.getBranch();
        final boolean branchDisp = !hasNaN(branch.getRawDisplacement().getData());
        final boolean branchAct = !hasNaN(branch.getRawActualization().getData());
        final boolean branchPrec = !hasNaN(branch.getRawPrecision().getData());
        final boolean aboveMean = !hasNaN(statistics.getAbove().getRawMean().getData());
        final boolean abovePrec = !hasNaN(statistics.getAbove().getRawPrecision().getData());
        final boolean belowMean = !hasNaN(statistics.getBelow().getRawMean().getData());
        final boolean belowPrec = !hasNaN(statistics.getBelow().getRawPrecision().getData());
        if (!branchDisp || !branchAct || !branchPrec || !aboveMean || !abovePrec || !belowMean || !belowPrec) {
            throw new IllegalStateException(
                    "External branch statistics already contain NaN"
                            + " branchDisp=" + branchDisp
                            + " branchAct=" + branchAct
                            + " branchPrec=" + branchPrec
                            + " aboveMean=" + aboveMean
                            + " abovePrec=" + abovePrec
                            + " belowMean=" + belowMean
                            + " belowPrec=" + belowPrec);
        }
    }

    private void ensureParameterScratchCapacity(final int parameterDimension) {
        if (parameterRestoreScratch.length < parameterDimension) {
            throw new IllegalStateException("parameter scratch too small for requested parameter dimension");
        }
    }

    private static double multiplyEntry(final DenseMatrix64F left,
                                        final DenseMatrix64F right,
                                        final int row,
                                        final int col) {
        double sum = 0.0;
        final int inner = left.numCols;
        for (int k = 0; k < inner; ++k) {
            sum += left.unsafe_get(row, k) * right.unsafe_get(k, col);
        }
        return sum;
    }

    private static double multiplyEntryTranspose(final DenseMatrix64F left,
                                                 final DenseMatrix64F middle,
                                                 final DenseMatrix64F right,
                                                 final int row,
                                                 final int col) {
        double sum = 0.0;
        final int inner = left.numRows;
        for (int k = 0; k < inner; ++k) {
            final double leftValue = left.unsafe_get(k, row);
            for (int l = 0; l < inner; ++l) {
                sum += leftValue * middle.unsafe_get(k, l) * right.unsafe_get(l, col);
            }
        }
        return sum;
    }

    private static void zero(final double[] vector) {
        for (int i = 0; i < vector.length; ++i) {
            vector[i] = 0.0;
        }
    }

    private static void zero(final double[][] matrix) {
        for (int i = 0; i < matrix.length; ++i) {
            for (int j = 0; j < matrix[i].length; ++j) {
                matrix[i][j] = 0.0;
            }
        }
    }

    private static void transposeInto(final double[][] in,
                                      final double[][] out) {
        for (int i = 0; i < in.length; ++i) {
            for (int j = 0; j < in[i].length; ++j) {
                out[i][j] = in[j][i];
            }
        }
    }

    private static void copyInto(final double[][] in,
                                 final double[][] out) {
        for (int i = 0; i < in.length; ++i) {
            System.arraycopy(in[i], 0, out[i], 0, in[i].length);
        }
    }

    private static void copyFromFlatInto(final double[] in, final double[][] out, final int dim) {
        for (int i = 0; i < dim; ++i) {
            System.arraycopy(in, i * dim, out[i], 0, dim);
        }
    }

    private static void transposeFromFlatInto(final double[] in, final double[][] out, final int dim) {
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                out[j][i] = in[i * dim + j];
            }
        }
    }

    private static void transposeFlatSquare(final double[] in, final double[] out, final int dim) {
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                out[j * dim + i] = in[i * dim + j];
            }
        }
    }

    private static void copyAdjoints(final CanonicalLocalTransitionAdjoints source,
                                     final CanonicalLocalTransitionAdjoints target,
                                     final boolean includeTransition,
                                     final boolean includeCovariance) {
        final int matLen = source.dLogL_dF.length;
        final int vecLen = source.dLogL_df.length;
        for (int k = 0; k < matLen; ++k) {
            target.dLogL_dF[k] = includeTransition ? source.dLogL_dF[k] : 0.0;
            target.dLogL_dOmega[k] = includeCovariance ? source.dLogL_dOmega[k] : 0.0;
        }
        for (int i = 0; i < vecLen; ++i) {
            target.dLogL_df[i] = includeTransition ? source.dLogL_df[i] : 0.0;
        }
    }

    private boolean useNumericLocalSelectionGradient() {
        return fallbackPolicy.useNumericLocalSelectionGradient();
    }

    private void emitDenseLocalSelectionDebug(final double branchLength,
                                              final double[] optimum,
                                              final CanonicalLocalTransitionAdjoints adjoints,
                                              final double[] analytic) {
        final double[] numeric = numericalLocalSelectionGradientDenseA(branchLength, optimum, adjoints);
        double maxAbs = 0.0;
        int maxIdx = -1;
        for (int i = 0; i < Math.min(analytic.length, numeric.length); ++i) {
            final double diff = Math.abs(analytic[i] - numeric[i]);
            if (diff > maxAbs) {
                maxAbs = diff;
                maxIdx = i;
            }
        }
        System.err.println("branchLocalSelectionFDDebug len=" + branchLength
                + " maxAbsDiff=" + maxAbs
                + " maxIdx=" + maxIdx
                + " analytic=" + Arrays.toString(analytic)
                + " numeric=" + Arrays.toString(numeric));
    }

    private void maybeEmitLocalMeanFDDebug(final double branchLength,
                                           final double[] optimum,
                                           final BranchSufficientStatistics statistics,
                                           final int nodeNumber,
                                           final double[] analytic) {
        if (!debugOptions.isBranchLocalMeanFiniteDifferenceEnabled()) {
            return;
        }
        final Integer requestedNode = debugOptions.getBranchLocalMeanFiniteDifferenceNode();
        if (requestedNode != null) {
            if (nodeNumber != requestedNode) {
                return;
            }
        }
        final double[] numericAdjointScore = numericalFrozenLocalMeanGradient(branchLength, optimum);
        final double[] numericFrozenFactor = numericalFrozenLocalMeanGradientFromFrozenFactor(
                branchLength, optimum, statistics);
        final double baseAdjointScore = evaluateLocalSelectionScore(branchLength, optimum, localAdjoints);
        final double baseFrozenFactor = branchWiring.evaluateFrozenLocalLogFactor(
                branchLength,
                optimum,
                statistics.getAbove(),
                statistics.getAboveParent(),
                statistics.getBelow());
        double maxAbsAdjoint = 0.0;
        int maxIdxAdjoint = -1;
        double maxAbsFrozenFactor = 0.0;
        int maxIdxFrozenFactor = -1;
        for (int i = 0; i < analytic.length; ++i) {
            final double adjointDiff = Math.abs(analytic[i] - numericAdjointScore[i]);
            if (adjointDiff > maxAbsAdjoint) {
                maxAbsAdjoint = adjointDiff;
                maxIdxAdjoint = i;
            }
            final double frozenFactorDiff = Math.abs(analytic[i] - numericFrozenFactor[i]);
            if (frozenFactorDiff > maxAbsFrozenFactor) {
                maxAbsFrozenFactor = frozenFactorDiff;
                maxIdxFrozenFactor = i;
            }
        }
        System.err.println("branchLocalMeanFDDebug node=" + nodeNumber
                + " len=" + branchLength
                + " belowPrecisionHasInf=" + hasInfinity(statsBelowPrecision(statistics))
                + " aboveParentPrecisionHasInf=" + hasInfinity(statsAboveParentPrecision(statistics))
                + " belowHasMissingMask=" + hasMissingMask(statistics)
                + " maxAbsAdjointScoreDiff=" + maxAbsAdjoint
                + " maxIdxAdjointScore=" + maxIdxAdjoint
                + " maxAbsFrozenFactorDiff=" + maxAbsFrozenFactor
                + " maxIdxFrozenFactor=" + maxIdxFrozenFactor
                + " baseScoreDelta=" + (baseAdjointScore - baseFrozenFactor)
                + " analytic=" + Arrays.toString(analytic)
                + " numericAdjointScore=" + Arrays.toString(numericAdjointScore)
                + " numericFrozenFactor=" + Arrays.toString(numericFrozenFactor));
    }

    private double[] numericalFrozenLocalMeanGradient(final double branchLength,
                                                      final double[] optimum) {
        final double[] gradient = new double[dimension];
        final double[] optimumWork = optimum.clone();
        for (int i = 0; i < dimension; ++i) {
            final double base = optimumWork[i];
            final double step = 1.0e-6 * Math.max(1.0, Math.abs(base));
            optimumWork[i] = base + step;
            final double plus = evaluateLocalSelectionScore(branchLength, optimumWork, localAdjoints);
            optimumWork[i] = base - step;
            final double minus = evaluateLocalSelectionScore(branchLength, optimumWork, localAdjoints);
            gradient[i] = (plus - minus) / (2.0 * step);
            optimumWork[i] = base;
        }
        return gradient;
    }

    private double[] numericalFrozenLocalMeanGradientFromFrozenFactor(final double branchLength,
                                                                      final double[] optimum,
                                                                      final BranchSufficientStatistics statistics) {
        final double[] gradient = new double[dimension];
        final double[] optimumWork = optimum.clone();
        for (int i = 0; i < dimension; ++i) {
            final double base = optimumWork[i];
            final double step = 1.0e-6 * Math.max(1.0, Math.abs(base));
            optimumWork[i] = base + step;
            final double plus = branchWiring.evaluateFrozenLocalLogFactor(
                    branchLength,
                    optimumWork,
                    statistics.getAbove(),
                    statistics.getAboveParent(),
                    statistics.getBelow());
            optimumWork[i] = base - step;
            final double minus = branchWiring.evaluateFrozenLocalLogFactor(
                    branchLength,
                    optimumWork,
                    statistics.getAbove(),
                    statistics.getAboveParent(),
                    statistics.getBelow());
            gradient[i] = (plus - minus) / (2.0 * step);
            optimumWork[i] = base;
        }
        return gradient;
    }

    private static DenseMatrix64F statsBelowPrecision(final BranchSufficientStatistics statistics) {
        return statistics.getBelow() == null ? null : statistics.getBelow().getRawPrecision();
    }

    private static DenseMatrix64F statsAboveParentPrecision(final BranchSufficientStatistics statistics) {
        return statistics.getAboveParent() == null ? null : statistics.getAboveParent().getRawPrecision();
    }

    private static boolean hasMissingMask(final BranchSufficientStatistics statistics) {
        final int[] missing = statistics.getMissing();
        if (missing == null) {
            return false;
        }
        for (int value : missing) {
            if (value != 0) {
                return true;
            }
        }
        return false;
    }

    private static boolean hasInfinity(final DenseMatrix64F matrix) {
        if (matrix == null) {
            return false;
        }
        final double[] data = matrix.getData();
        for (double value : data) {
            if (Double.isInfinite(value)) {
                return true;
            }
        }
        return false;
    }

    private double[] numericalLocalSelectionGradientDenseA(final double branchLength,
                                                           final double[] optimum,
                                                           final CanonicalLocalTransitionAdjoints adjoints) {
        final int d = dimension;
        final double[] gradient = new double[d * d];
        final double[] restore = new double[d * d];
        final dr.inference.model.MatrixParameterInterface drift = processModel.getDriftMatrix();

        int k = 0;
        for (int row = 0; row < d; ++row) {
            for (int col = 0; col < d; ++col) {
                restore[k++] = drift.getParameterValue(row, col);
            }
        }

        try {
            k = 0;
            for (int row = 0; row < d; ++row) {
                for (int col = 0; col < d; ++col) {
                    final double base = restore[k];
                    final double step = 1.0e-6 * Math.max(1.0, Math.abs(base));
                    drift.setParameterValue(row, col, base + step);
                    final double plus = evaluateLocalSelectionScore(branchLength, optimum, adjoints);
                    drift.setParameterValue(row, col, base - step);
                    final double minus = evaluateLocalSelectionScore(branchLength, optimum, adjoints);
                    gradient[k] = (plus - minus) / (2.0 * step);
                    drift.setParameterValue(row, col, base);
                    ++k;
                }
            }
            return gradient;
        } finally {
            k = 0;
            for (int row = 0; row < d; ++row) {
                for (int col = 0; col < d; ++col) {
                    final double original = restore[k++];
                    if (drift.getParameterValue(row, col) != original) {
                        drift.setParameterValue(row, col, original);
                    }
                }
            }
        }
    }

    private double[] numericalLocalSelectionGradientFromFrozenFactor(final double branchLength,
                                                                     final double[] optimum,
                                                                     final BranchSufficientStatistics statistics) {
        // Keep FD objective aligned with analytic selection gradient:
        // finite-difference the same local adjoint score used by the analytic chain rule.
        branchWiring.fillLocalAdjoints(branchLength, optimum, statistics, localAdjoints);
        return numericalLocalSelectionGradientDenseA(branchLength, optimum, localAdjoints);
    }

    private double[] assembleBlockGradientResult(final Parameter requestedParameter,
                                                 final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter,
                                                 final double[] nativeBlockGradient,
                                                 final double[][] gradientR) {
        if (requestedParameter == blockParameter.getParameter()) {
            if (blockParameter.getRotationMatrixParameter() instanceof OrthogonalMatrixProvider) {
                final double[] angleGradient =
                        ((OrthogonalMatrixProvider) blockParameter.getRotationMatrixParameter()).pullBackGradient(gradientR);
                final double[] out = new double[blockParameter.getBlockDiagonalNParameters() + angleGradient.length];
                System.arraycopy(nativeBlockGradient, 0, out, 0, blockParameter.getBlockDiagonalNParameters());
                System.arraycopy(angleGradient, 0, out, blockParameter.getBlockDiagonalNParameters(), angleGradient.length);
                return out;
            }
            final double[] out = new double[blockParameter.getBlockDiagonalNParameters() + dimension * dimension];
            System.arraycopy(nativeBlockGradient, 0, out, 0, blockParameter.getBlockDiagonalNParameters());
            final double[] flatRotation = flattenRotationGradient(blockParameter.getRotationMatrixParameter(), gradientR);
            System.arraycopy(flatRotation, 0, out, blockParameter.getBlockDiagonalNParameters(), flatRotation.length);
            return out;
        }

        if (requestedParameter == blockParameter.getRotationMatrixParameter()) {
            return flattenRotationGradient(blockParameter.getRotationMatrixParameter(), gradientR);
        }
        if (blockParameter.getRotationMatrixParameter() instanceof OrthogonalMatrixProvider
                && requestedParameter == ((OrthogonalMatrixProvider) blockParameter.getRotationMatrixParameter()).getOrthogonalParameter()) {
            return ((OrthogonalMatrixProvider) blockParameter.getRotationMatrixParameter()).pullBackGradient(gradientR);
        }
        if (requestedParameter == blockParameter.getScalarBlockParameter()) {
            return new double[]{nativeBlockGradient[0]};
        }

        final int blockBase = blockParameter.hasLeadingOneByOneBlock() ? 1 : 0;
        final int blockWidth = blockParameter.getNum2x2Blocks();
        for (int family = 0; family < blockParameter.getTwoByTwoParameterFamilyCount(); ++family) {
            if (requestedParameter == blockParameter.getTwoByTwoBlockParameter(family)) {
                final double[] out = new double[blockWidth];
                System.arraycopy(nativeBlockGradient, blockBase + family * blockWidth, out, 0, blockWidth);
                return out;
            }
        }
        throw new IllegalArgumentException("Unsupported block parameter: " + requestedParameter.getId());
    }

    private static double[] flattenColumnMajor(final double[][] matrix) {
        final int d = matrix.length;
        final double[] out = new double[d * d];
        int index = 0;
        for (int col = 0; col < d; ++col) {
            for (int row = 0; row < d; ++row) {
                out[index++] = matrix[row][col];
            }
        }
        return out;
    }

    private static double[] flattenRowMajor(final double[][] matrix) {
        final int d = matrix.length;
        final double[] out = new double[d * d];
        int index = 0;
        for (int row = 0; row < d; ++row) {
            for (int col = 0; col < d; ++col) {
                out[index++] = matrix[row][col];
            }
        }
        return out;
    }

    private static double[] flattenRotationGradient(final Parameter rotationParameter,
                                                    final double[][] gradientR) {
        if (rotationParameter instanceof TransposedMatrixParameter) {
            return flattenRowMajor(gradientR);
        }
        return flattenColumnMajor(gradientR);
    }

    private double[] reorderDenseGradientForRequestedParameter(final Parameter requestedParameter,
                                                               final double[] rowMajorGradient) {
        if (requestedParameter instanceof TransposedMatrixParameter) {
            return rowMajorGradient;
        }
        if (!(requestedParameter instanceof dr.inference.model.MatrixParameterInterface)) {
            return rowMajorGradient;
        }
        final double[] reordered = new double[rowMajorGradient.length];
        int index = 0;
        for (int col = 0; col < dimension; ++col) {
            for (int row = 0; row < dimension; ++row) {
                reordered[index++] = rowMajorGradient[row * dimension + col];
            }
        }
        return reordered;
    }

    private double[] pullBackDenseGradientToBlock(final Parameter requestedParameter,
                                                  final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter,
                                                  final double[] denseGradient) {
        final int d = dimension;
        final double[][] gradientA = new double[d][d];
        final double[] rData = new double[d * d];
        final double[] rinvData = new double[d * d];
        final double[] rawBlockSource = new double[blockParameter.getCompressedDDimension()];
        final double[] nativeBlockGradient = new double[blockParameter.getBlockDiagonalNParameters()];
        final double[][] gradientD = new double[d][d];
        final double[][] gradientR = new double[d][d];
        final double[][] aMatrix = blockParameter.getParameterAsMatrix();

        fillMatrixFromRowMajorGradient(gradientA, denseGradient, d);
        blockParameter.fillRAndRinv(rData, rinvData);
        computeGradientWrtBlockDiagonalBasis(gradientA, rData, rinvData, gradientD);
        computeGradientWrtRotationMatrix(gradientA, aMatrix, rData, rinvData, gradientR);
        compressActiveBlockGradient(blockParameter, gradientD, rawBlockSource);
        blockParameter.chainGradient(rawBlockSource, nativeBlockGradient);
        return assembleBlockGradientResult(requestedParameter, blockParameter, nativeBlockGradient, gradientR);
    }

    private static void fillMatrixFromRowMajorGradient(final double[][] out,
                                                       final double[] gradient,
                                                       final int d) {
        for (int row = 0; row < d; ++row) {
            for (int col = 0; col < d; ++col) {
                out[row][col] = gradient[row * d + col];
            }
        }
    }

    private static void computeGradientWrtBlockDiagonalBasis(final double[][] gradientA,
                                                             final double[] rData,
                                                             final double[] rinvData,
                                                             final double[][] out) {
        final int d = gradientA.length;
        final double[][] temp = new double[d][d];
        final double[][] rinvTranspose = new double[d][d];

        for (int row = 0; row < d; ++row) {
            for (int col = 0; col < d; ++col) {
                rinvTranspose[row][col] = rinvData[col * d + row];
            }
        }

        multiplyTransposeLeft(rData, d, gradientA, temp);
        multiplyDense(temp, rinvTranspose, out);
    }

    private static void computeGradientWrtRotationMatrix(final double[][] gradientA,
                                                         final double[][] aMatrix,
                                                         final double[] rData,
                                                         final double[] rinvData,
                                                         final double[][] out) {
        final int d = gradientA.length;
        final double[][] aTranspose = new double[d][d];
        final double[][] dMatrix = new double[d][d];
        final double[][] dTranspose = new double[d][d];
        final double[][] temp = new double[d][d];
        final double[][] temp2 = new double[d][d];
        final double[][] rinv = new double[d][d];
        final double[][] rinvTranspose = new double[d][d];
        final double[][] rMatrix = new double[d][d];

        transpose(aMatrix, aTranspose);

        for (int row = 0; row < d; ++row) {
            for (int col = 0; col < d; ++col) {
                rinv[row][col] = rinvData[row * d + col];
                rinvTranspose[row][col] = rinvData[col * d + row];
                rMatrix[row][col] = rData[row * d + col];
            }
        }

        multiplyDense(rinv, aMatrix, temp2);
        multiplyDense(temp2, rMatrix, dMatrix);
        transpose(dMatrix, dTranspose);

        multiplyDense(gradientA, rinvTranspose, temp);
        multiplyDense(temp, dTranspose, temp2);

        multiplyDense(aTranspose, gradientA, temp);
        multiplyDense(temp, rinvTranspose, out);

        subtractIntoOut(temp2, out);
    }

    private static void compressActiveBlockGradient(final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter,
                                                    final double[][] gradientD,
                                                    final double[] out) {
        final int d = gradientD.length;
        int blockIndex = 0;
        final int upperBase = d;
        final int lowerBase = d + blockParameter.getNum2x2Blocks();

        for (int i = 0; i < d; ++i) {
            out[i] = gradientD[i][i];
        }

        for (int b = 0; b < blockParameter.getNumBlocks(); ++b) {
            if (blockParameter.getBlockSizes()[b] != 2) {
                continue;
            }
            final int start = blockParameter.getBlockStarts()[b];
            out[upperBase + blockIndex] = gradientD[start][start + 1];
            out[lowerBase + blockIndex] = gradientD[start + 1][start];
            blockIndex++;
        }
    }

    private static void multiplyTransposeLeft(final double[] leftRowMajor,
                                              final int dimension,
                                              final double[][] right,
                                              final double[][] out) {
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                double sum = 0.0;
                for (int k = 0; k < dimension; ++k) {
                    sum += leftRowMajor[k * dimension + i] * right[k][j];
                }
                out[i][j] = sum;
            }
        }
    }

    private static void multiplyDense(final double[][] left,
                                      final double[][] right,
                                      final double[][] out) {
        final int rows = left.length;
        final int cols = right[0].length;
        final int inner = right.length;
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                double sum = 0.0;
                for (int k = 0; k < inner; ++k) {
                    sum += left[i][k] * right[k][j];
                }
                out[i][j] = sum;
            }
        }
    }

    private static void transpose(final double[][] in,
                                  final double[][] out) {
        for (int i = 0; i < in.length; ++i) {
            for (int j = 0; j < in[i].length; ++j) {
                out[j][i] = in[i][j];
            }
        }
    }

    private static void subtractIntoOut(final double[][] left,
                                        final double[][] out) {
        for (int i = 0; i < left.length; ++i) {
            for (int j = 0; j < left[i].length; ++j) {
                out[i][j] = left[i][j] - out[i][j];
            }
        }
    }

}
