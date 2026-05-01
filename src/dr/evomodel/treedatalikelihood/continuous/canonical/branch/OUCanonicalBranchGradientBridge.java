package dr.evomodel.treedatalikelihood.continuous.canonical.branch;

import dr.evomodel.treedatalikelihood.continuous.CanonicalDebugOptions;
import dr.evomodel.treedatalikelihood.continuous.CanonicalGradientFallbackPolicy;
import dr.evomodel.treedatalikelihood.continuous.OUGaussianBranchTransitionProvider;

import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalLocalTransitionAdjoints;
import dr.evomodel.treedatalikelihood.preorder.BranchSufficientStatistics;
import dr.inference.model.AbstractBlockDiagonalTwoByTwoMatrixParameter;
import dr.inference.model.Parameter;

/**
 * Tree-side facade for canonical OU branch gradient helpers.
 */
public final class OUCanonicalBranchGradientBridge {
    private static final ThreadLocal<Integer> DEBUG_NODE_CONTEXT = new ThreadLocal<Integer>();

    private final OUGaussianBranchTransitionProvider branchTransitionProvider;
    private final OUCanonicalBranchWiring branchWiring;
    private final int dimension;
    private final OUCanonicalBranchGradientScratch scratch;
    private final OUCanonicalBranchMeanGradient meanGradient;
    private final OUCanonicalBranchDiffusionGradient diffusionGradient;
    private final OUCanonicalBranchSelectionGradient selectionGradient;

    public OUCanonicalBranchGradientBridge(
            final OUGaussianBranchTransitionProvider branchTransitionProvider) {
        this(branchTransitionProvider,
                CanonicalDebugOptions.fromSystemProperties(),
                CanonicalGradientFallbackPolicy.fromSystemProperties());
    }

    OUCanonicalBranchGradientBridge(
            final OUGaussianBranchTransitionProvider branchTransitionProvider,
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
        this.branchWiring = new OUCanonicalBranchWiring(branchTransitionProvider, debugOptions, fallbackPolicy);
        final OUProcessModel processModel = branchWiring.getProcessModel();
        this.dimension = branchWiring.getDimension();
        this.scratch = new OUCanonicalBranchGradientScratch(dimension);
        final OUCanonicalBranchFiniteDifference finiteDifference =
                new OUCanonicalBranchFiniteDifference(
                        branchTransitionProvider,
                        branchWiring,
                        processModel,
                        debugOptions,
                        dimension,
                        scratch);
        final OUCanonicalNativeSelectionGradient nativeSelectionGradient =
                new OUCanonicalNativeSelectionGradient(
                        branchWiring,
                        processModel,
                        debugOptions,
                        fallbackPolicy,
                        finiteDifference,
                        scratch,
                        dimension,
                        DEBUG_NODE_CONTEXT);
        this.meanGradient = new OUCanonicalBranchMeanGradient(branchWiring, finiteDifference, scratch);
        this.diffusionGradient = new OUCanonicalBranchDiffusionGradient(
                branchTransitionProvider,
                branchWiring,
                processModel,
                dimension,
                scratch);
        this.selectionGradient = new OUCanonicalBranchSelectionGradient(
                branchTransitionProvider,
                branchWiring,
                processModel,
                debugOptions,
                fallbackPolicy,
                finiteDifference,
                nativeSelectionGradient,
                scratch,
                dimension);
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
        refreshProcessSnapshots();
        return meanGradient.getGradientWrtBranchDrift(statistics);
    }

    public double[] getGradientWrtBranchDrift(final double branchLength,
                                              final double[] optimum,
                                              final BranchSufficientStatistics statistics,
                                              final int nodeNumber) {
        refreshProcessSnapshots();
        return meanGradient.getGradientWrtBranchDrift(branchLength, optimum, statistics, nodeNumber);
    }

    public double[] getGradientWrtVariance(final double branchLength,
                                           final double[] optimum,
                                           final BranchSufficientStatistics statistics) {
        return diffusionGradient.getGradientWrtVariance(branchLength, optimum, statistics);
    }

    public void accumulateGlobalRemainderSelectionGradientForBranch(
            final double branchLength,
            final double[] barVdiFlatIn,
            final AbstractBlockDiagonalTwoByTwoMatrixParameter nativeBlockParameter,
            final Parameter requestedParameter,
            final double[] gradientAccumulator) {
        selectionGradient.accumulateGlobalRemainderSelectionGradientForBranch(
                branchLength,
                barVdiFlatIn,
                nativeBlockParameter,
                requestedParameter,
                gradientAccumulator);
    }

    public double[] getGradientWrtSelection(final double branchLength,
                                            final double[] optimum,
                                            final BranchSufficientStatistics statistics) {
        return selectionGradient.getGradientWrtSelection(branchLength, optimum, statistics);
    }

    public double[] getGradientWrtSelectionComponents(final double branchLength,
                                                      final double[] optimum,
                                                      final BranchSufficientStatistics statistics) {
        return selectionGradient.getGradientWrtSelectionComponents(branchLength, optimum, statistics);
    }

    public double[] getNumericalDenseGradientWrtSelectionComponents(final double branchLength,
                                                                    final double[] optimum,
                                                                    final BranchSufficientStatistics statistics) {
        return selectionGradient.getNumericalDenseGradientWrtSelectionComponents(branchLength, optimum, statistics);
    }

    public double[] getGradientWrtSelection(final double branchLength,
                                            final double[] optimum,
                                            final BranchSufficientStatistics statistics,
                                            final boolean externalBranch,
                                            final AbstractBlockDiagonalTwoByTwoMatrixParameter nativeBlockParameter,
                                            final Parameter requestedParameter) {
        return selectionGradient.getGradientWrtSelection(
                branchLength,
                optimum,
                statistics,
                externalBranch,
                nativeBlockParameter,
                requestedParameter);
    }

    public double[] getNumericalFrozenLocalGradientWrtSelection(
            final double branchLength,
            final double[] optimum,
            final BranchSufficientStatistics statistics,
            final AbstractBlockDiagonalTwoByTwoMatrixParameter nativeBlockParameter,
            final Parameter requestedParameter) {
        return selectionGradient.getNumericalFrozenLocalGradientWrtSelection(
                branchLength,
                optimum,
                statistics,
                nativeBlockParameter,
                requestedParameter);
    }

    private void refreshProcessSnapshots() {
        branchTransitionProvider.getProcessModel();
    }
}
