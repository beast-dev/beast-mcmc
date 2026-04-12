package dr.evomodel.treedatalikelihood.continuous.cdi;

import dr.evolution.tree.NodeRef;
import dr.evomodel.continuous.MultivariateElasticModel;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.continuous.OULegacyAttenuationGradientHelper;
import dr.evomodel.treedatalikelihood.continuous.OULegacyActualizationGradientHelper;
import dr.evomodel.treedatalikelihood.continuous.OUDiffusionModelDelegate;
import dr.evomodel.treedatalikelihood.preorder.BranchSufficientStatistics;
import org.ejml.data.DenseMatrix64F;
import java.util.Arrays;

public final class OUActualizationStrategies {

    private static final OUGradientStrategy DIAGONAL_GRADIENT = new OUGradientStrategy() {
        @Override
        public DenseMatrix64F getGradientVarianceWrtVariance(final OUDiffusionModelDelegate delegate,
                                                             final NodeRef node,
                                                             final ContinuousDiffusionIntegrator cdi,
                                                             final ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                                             final DenseMatrix64F gradient) {
            if (delegate.getTree().isRoot(node)) {
                return delegate.rootGradientVarianceWrtVariance(node, cdi, likelihoodDelegate, gradient);
            }
            final DenseMatrix64F result = gradient.copy();
            OULegacyActualizationGradientHelper.actualizeGradientDiagonal(
                    delegate.getElasticModel(),
                    delegate.getDimTrait(),
                    cdi,
                    delegate.getMatrixBufferOffsetIndex(node.getNumber()),
                    node.getNumber(),
                    result);
            return result;
        }

        @Override
        public DenseMatrix64F getGradientVarianceWrtAttenuation(final OUDiffusionModelDelegate delegate,
                                                                final NodeRef node,
                                                                final ContinuousDiffusionIntegrator cdi,
                                                                final BranchSufficientStatistics statistics,
                                                                final DenseMatrix64F gradient) {
            if (delegate.getTree().isRoot(node)) {
                throw new AssertionError("Gradient wrt actualization is not available for the root.");
            }
            return OULegacyAttenuationGradientHelper.getVarianceWrtAttenuationDiagonal(
                    delegate, node, cdi, statistics, gradient);
        }

        @Override
        public DenseMatrix64F getGradientDisplacementWrtAttenuation(final OUDiffusionModelDelegate delegate,
                                                                    final NodeRef node,
                                                                    final ContinuousDiffusionIntegrator cdi,
                                                                    final BranchSufficientStatistics statistics,
                                                                    final DenseMatrix64F gradient) {
            if (delegate.getTree().isRoot(node)) {
                throw new AssertionError("Gradient wrt actualization is not available for the root.");
            }
            return OULegacyAttenuationGradientHelper.getDisplacementWrtAttenuationDiagonal(
                    delegate, node, cdi, statistics, gradient);
        }

        @Override
        public DenseMatrix64F getGradientDisplacementWrtDrift(final OUDiffusionModelDelegate delegate,
                                                              final NodeRef node,
                                                              final ContinuousDiffusionIntegrator cdi,
                                                              final ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                                              final DenseMatrix64F gradient) {
            final DenseMatrix64F inputSnapshot = gradient.copy();
            final DenseMatrix64F result = gradient.copy();
            final int matrixBufferOffsetIndex = delegate.getMatrixBufferOffsetIndex(node.getNumber());
            if (delegate.usesDiagonalActualizationBuffer(cdi)) {
                OULegacyActualizationGradientHelper.actualizeDisplacementGradientDiagonal(
                        delegate.getDimTrait(),
                        cdi,
                        matrixBufferOffsetIndex,
                        result);
            } else {
                OULegacyActualizationGradientHelper.actualizeDisplacementGradient(
                        delegate.getDimTrait(),
                        cdi,
                        matrixBufferOffsetIndex,
                        result);
            }
            maybeEmitDriftActualizationDebug(delegate, node, cdi, matrixBufferOffsetIndex, inputSnapshot, result);
            return result;
        }

        @Override
        public double[] getGradientDisplacementWrtRoot(final OUDiffusionModelDelegate delegate,
                                                       final NodeRef node,
                                                       final ContinuousDiffusionIntegrator cdi,
                                                       final ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                                       final DenseMatrix64F gradient) {
            final boolean fixedRoot = delegate.hasFixedRoot(likelihoodDelegate);
            if (fixedRoot && delegate.getTree().isRoot(delegate.getTree().getParent(node))) {
                if (delegate.usesDiagonalActualizationBuffer(cdi)) {
                    return OULegacyActualizationGradientHelper.actualizeRootGradientDiagonal(
                            delegate.getDimTrait(),
                            cdi,
                            delegate.getMatrixBufferOffsetIndex(node.getNumber()),
                            gradient);
                }
                return OULegacyActualizationGradientHelper.actualizeRootGradientFull(
                        delegate.getDimTrait(),
                        cdi,
                        delegate.getMatrixBufferOffsetIndex(node.getNumber()),
                        gradient);
            }
            if (!fixedRoot && delegate.getTree().isRoot(node)) {
                return gradient.getData();
            }
            return new double[gradient.getNumRows()];
        }
    };

    private static final OUGradientStrategy GENERAL_GRADIENT = new OUGradientStrategy() {
        @Override
        public DenseMatrix64F getGradientVarianceWrtVariance(final OUDiffusionModelDelegate delegate,
                                                             final NodeRef node,
                                                             final ContinuousDiffusionIntegrator cdi,
                                                             final ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                                             final DenseMatrix64F gradient) {
            if (delegate.getTree().isRoot(node)) {
                return delegate.rootGradientVarianceWrtVariance(node, cdi, likelihoodDelegate, gradient);
            }
            final DenseMatrix64F result = gradient.copy();
            OULegacyActualizationGradientHelper.actualizeGradient(
                    delegate.getElasticModel(),
                    delegate.getDimTrait(),
                    cdi,
                    delegate.getMatrixBufferOffsetIndex(node.getNumber()),
                    node.getNumber(),
                    result);
            return result;
        }

        @Override
        public DenseMatrix64F getGradientVarianceWrtAttenuation(final OUDiffusionModelDelegate delegate,
                                                                final NodeRef node,
                                                                final ContinuousDiffusionIntegrator cdi,
                                                                final BranchSufficientStatistics statistics,
                                                                final DenseMatrix64F gradient) {
            return OULegacyAttenuationGradientHelper.unsupportedVarianceWrtAttenuation("general");
        }

        @Override
        public DenseMatrix64F getGradientDisplacementWrtAttenuation(final OUDiffusionModelDelegate delegate,
                                                                    final NodeRef node,
                                                                    final ContinuousDiffusionIntegrator cdi,
                                                                    final BranchSufficientStatistics statistics,
                                                                    final DenseMatrix64F gradient) {
            return OULegacyAttenuationGradientHelper.unsupportedDisplacementWrtAttenuation("general");
        }

        @Override
        public DenseMatrix64F getGradientDisplacementWrtDrift(final OUDiffusionModelDelegate delegate,
                                                              final NodeRef node,
                                                              final ContinuousDiffusionIntegrator cdi,
                                                              final ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                                              final DenseMatrix64F gradient) {
            final DenseMatrix64F inputSnapshot = gradient.copy();
            final DenseMatrix64F result = gradient.copy();
            final int matrixBufferOffsetIndex = delegate.getMatrixBufferOffsetIndex(node.getNumber());
            OULegacyActualizationGradientHelper.actualizeDisplacementGradient(
                    delegate.getDimTrait(),
                    cdi,
                    matrixBufferOffsetIndex,
                    result);
            maybeEmitDriftActualizationDebug(delegate, node, cdi, matrixBufferOffsetIndex, inputSnapshot, result);
            return result;
        }

        @Override
        public double[] getGradientDisplacementWrtRoot(final OUDiffusionModelDelegate delegate,
                                                       final NodeRef node,
                                                       final ContinuousDiffusionIntegrator cdi,
                                                       final ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                                       final DenseMatrix64F gradient) {
            final boolean fixedRoot = delegate.hasFixedRoot(likelihoodDelegate);
            if (fixedRoot && delegate.getTree().isRoot(delegate.getTree().getParent(node))) {
                return OULegacyActualizationGradientHelper.actualizeRootGradientFull(
                        delegate.getDimTrait(),
                        cdi,
                        delegate.getMatrixBufferOffsetIndex(node.getNumber()),
                        gradient);
            }
            if (!fixedRoot && delegate.getTree().isRoot(node)) {
                return gradient.getData();
            }
            return new double[gradient.getNumRows()];
        }
    };

    private static final OUGradientStrategy BLOCK_GRADIENT = new OUGradientStrategy() {
        @Override
        public DenseMatrix64F getGradientVarianceWrtVariance(final OUDiffusionModelDelegate delegate,
                                                             final NodeRef node,
                                                             final ContinuousDiffusionIntegrator cdi,
                                                             final ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                                             final DenseMatrix64F gradient) {
            if (delegate.getTree().isRoot(node)) {
                return delegate.rootGradientVarianceWrtVariance(node, cdi, likelihoodDelegate, gradient);
            }
            final DenseMatrix64F result = gradient.copy();
            if (delegate.usesDiagonalActualizationBuffer(cdi)) {
                OULegacyActualizationGradientHelper.actualizeGradientDiagonal(
                        delegate.getElasticModel(),
                        delegate.getDimTrait(),
                        cdi,
                        delegate.getMatrixBufferOffsetIndex(node.getNumber()),
                        node.getNumber(),
                        result);
            } else {
                OULegacyActualizationGradientHelper.actualizeGradient(
                        delegate.getElasticModel(),
                        delegate.getDimTrait(),
                        cdi,
                        delegate.getMatrixBufferOffsetIndex(node.getNumber()),
                        node.getNumber(),
                        result);
            }
            return result;
        }

        @Override
        public DenseMatrix64F getGradientVarianceWrtAttenuation(final OUDiffusionModelDelegate delegate,
                                                                final NodeRef node,
                                                                final ContinuousDiffusionIntegrator cdi,
                                                                final BranchSufficientStatistics statistics,
                                                                final DenseMatrix64F gradient) {
            return OULegacyAttenuationGradientHelper.unsupportedVarianceWrtAttenuation("block");
        }

        @Override
        public DenseMatrix64F getGradientDisplacementWrtAttenuation(final OUDiffusionModelDelegate delegate,
                                                                    final NodeRef node,
                                                                    final ContinuousDiffusionIntegrator cdi,
                                                                    final BranchSufficientStatistics statistics,
                                                                    final DenseMatrix64F gradient) {
            return OULegacyAttenuationGradientHelper.unsupportedDisplacementWrtAttenuation("block");
        }

        @Override
        public DenseMatrix64F getGradientDisplacementWrtDrift(final OUDiffusionModelDelegate delegate,
                                                              final NodeRef node,
                                                              final ContinuousDiffusionIntegrator cdi,
                                                              final ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                                              final DenseMatrix64F gradient) {
            final DenseMatrix64F result = gradient.copy();
            if (delegate.usesDiagonalActualizationBuffer(cdi)) {
                OULegacyActualizationGradientHelper.actualizeDisplacementGradientDiagonal(
                        delegate.getDimTrait(),
                        cdi,
                        delegate.getMatrixBufferOffsetIndex(node.getNumber()),
                        result);
            } else {
                OULegacyActualizationGradientHelper.actualizeDisplacementGradient(
                        delegate.getDimTrait(),
                        cdi,
                        delegate.getMatrixBufferOffsetIndex(node.getNumber()),
                        result);
            }
            return result;
        }

        @Override
        public double[] getGradientDisplacementWrtRoot(final OUDiffusionModelDelegate delegate,
                                                       final NodeRef node,
                                                       final ContinuousDiffusionIntegrator cdi,
                                                       final ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                                       final DenseMatrix64F gradient) {
            final boolean fixedRoot = delegate.hasFixedRoot(likelihoodDelegate);
            if (fixedRoot && delegate.getTree().isRoot(delegate.getTree().getParent(node))) {
                if (delegate.usesDiagonalActualizationBuffer(cdi)) {
                    return OULegacyActualizationGradientHelper.actualizeRootGradientDiagonal(
                            delegate.getDimTrait(),
                            cdi,
                            delegate.getMatrixBufferOffsetIndex(node.getNumber()),
                            gradient);
                }
                return OULegacyActualizationGradientHelper.actualizeRootGradientFull(
                        delegate.getDimTrait(),
                        cdi,
                        delegate.getMatrixBufferOffsetIndex(node.getNumber()),
                        gradient);
            }
            if (!fixedRoot && delegate.getTree().isRoot(node)) {
                return gradient.getData();
            }
            return new double[gradient.getNumRows()];
        }
    };

    private static final OUActualizationStrategy DIAGONAL = new OUActualizationStrategy() {
        @Override
        public void setDiffusionStationaryVariance(final SafeMultivariateActualizedWithDriftIntegrator integrator,
                                                   final int precisionIndex,
                                                   final double[] basisValues,
                                                   final double[] basisRotations) {
            integrator.performDiagonalSetDiffusionStationaryVariance(precisionIndex, basisValues, basisRotations);
        }

        @Override
        public void updateOrnsteinUhlenbeckDiffusionMatrices(final SafeMultivariateActualizedWithDriftIntegrator integrator,
                                                             final int precisionIndex,
                                                             final int[] probabilityIndices,
                                                             final double[] edgeLengths,
                                                             final double[] optimalRates,
                                                             final double[] basisValues,
                                                             final double[] basisRotations,
                                                             final int updateCount) {
            integrator.performDiagonalUpdateOrnsteinUhlenbeckDiffusionMatrices(
                    precisionIndex, probabilityIndices, edgeLengths, optimalRates, basisValues, basisRotations, updateCount);
        }

        @Override
        public void updateIntegratedOrnsteinUhlenbeckDiffusionMatrices(final SafeMultivariateActualizedWithDriftIntegrator integrator,
                                                                       final int precisionIndex,
                                                                       final int[] probabilityIndices,
                                                                       final double[] edgeLengths,
                                                                       final double[] optimalRates,
                                                                       final double[] basisValues,
                                                                       final double[] basisRotations,
                                                                       final int updateCount) {
            integrator.performDiagonalUpdateIntegratedOrnsteinUhlenbeckDiffusionMatrices(
                    precisionIndex, probabilityIndices, edgeLengths, optimalRates, basisValues, basisRotations, updateCount);
        }
    };

    private static void maybeEmitDriftActualizationDebug(final OUDiffusionModelDelegate delegate,
                                                         final NodeRef node,
                                                         final ContinuousDiffusionIntegrator cdi,
                                                         final int matrixBufferOffsetIndex,
                                                         final DenseMatrix64F input,
                                                         final DenseMatrix64F output) {
        if (!Boolean.getBoolean("beast.debug.ouDriftActualization")) {
            return;
        }
        final String nodeFilter = System.getProperty("beast.debug.ouDriftActualization.node");
        if (nodeFilter != null && !nodeFilter.isEmpty()) {
            try {
                if (node.getNumber() != Integer.parseInt(nodeFilter.trim())) {
                    return;
                }
            } catch (NumberFormatException ignored) {
                // fall through and keep logging if the node filter is malformed
            }
        }

        final int dim = delegate.getDimTrait();
        final double[] branchOneMinusActualization = new double[dim * dim];
        cdi.getBranch1mActualization(matrixBufferOffsetIndex, branchOneMinusActualization);
        System.err.println("ouDriftActualization"
                + " node=" + node.getNumber()
                + " matrixBufferOffsetIndex=" + matrixBufferOffsetIndex
                + " in=" + Arrays.toString(input.getData())
                + " branch1mActualization=" + Arrays.toString(branchOneMinusActualization)
                + " out=" + Arrays.toString(output.getData()));
    }

    private static final OUActualizationStrategy GENERAL = new OUActualizationStrategy() {
        @Override
        public void setDiffusionStationaryVariance(final SafeMultivariateActualizedWithDriftIntegrator integrator,
                                                   final int precisionIndex,
                                                   final double[] basisValues,
                                                   final double[] basisRotations) {
            integrator.performGeneralSetDiffusionStationaryVariance(precisionIndex, basisValues, basisRotations);
        }

        @Override
        public void updateOrnsteinUhlenbeckDiffusionMatrices(final SafeMultivariateActualizedWithDriftIntegrator integrator,
                                                             final int precisionIndex,
                                                             final int[] probabilityIndices,
                                                             final double[] edgeLengths,
                                                             final double[] optimalRates,
                                                             final double[] basisValues,
                                                             final double[] basisRotations,
                                                             final int updateCount) {
            integrator.performGeneralUpdateOrnsteinUhlenbeckDiffusionMatrices(
                    precisionIndex, probabilityIndices, edgeLengths, optimalRates, basisValues, basisRotations, updateCount);
        }

        @Override
        public void updateIntegratedOrnsteinUhlenbeckDiffusionMatrices(final SafeMultivariateActualizedWithDriftIntegrator integrator,
                                                                       final int precisionIndex,
                                                                       final int[] probabilityIndices,
                                                                       final double[] edgeLengths,
                                                                       final double[] optimalRates,
                                                                       final double[] basisValues,
                                                                       final double[] basisRotations,
                                                                       final int updateCount) {
            integrator.performGeneralUpdateIntegratedOrnsteinUhlenbeckDiffusionMatrices(
                    precisionIndex, probabilityIndices, edgeLengths, optimalRates, basisValues, basisRotations, updateCount);
        }
    };

    private static final OUActualizationStrategy BLOCK = new OUActualizationStrategy() {
        @Override
        public void setDiffusionStationaryVariance(final SafeMultivariateActualizedWithDriftIntegrator integrator,
                                                   final int precisionIndex,
                                                   final double[] basisValues,
                                                   final double[] basisRotations) {
            integrator.performBlockSetDiffusionStationaryVariance(precisionIndex, basisValues, basisRotations);
        }

        @Override
        public void updateOrnsteinUhlenbeckDiffusionMatrices(final SafeMultivariateActualizedWithDriftIntegrator integrator,
                                                             final int precisionIndex,
                                                             final int[] probabilityIndices,
                                                             final double[] edgeLengths,
                                                             final double[] optimalRates,
                                                             final double[] basisValues,
                                                             final double[] basisRotations,
                                                             final int updateCount) {
            integrator.performBlockUpdateOrnsteinUhlenbeckDiffusionMatrices(
                    precisionIndex, probabilityIndices, edgeLengths, optimalRates, basisValues, basisRotations, updateCount);
        }

        @Override
        public void updateIntegratedOrnsteinUhlenbeckDiffusionMatrices(final SafeMultivariateActualizedWithDriftIntegrator integrator,
                                                                       final int precisionIndex,
                                                                       final int[] probabilityIndices,
                                                                       final double[] edgeLengths,
                                                                       final double[] optimalRates,
                                                                       final double[] basisValues,
                                                                       final double[] basisRotations,
                                                                       final int updateCount) {
            throw new UnsupportedOperationException(
                    "Integrated OU block basis is not yet implemented in SafeMultivariateActualizedWithDriftIntegrator");
        }
    };

    private OUActualizationStrategies() {
    }

    static OUActualizationStrategy general() {
        return GENERAL;
    }

    static OUActualizationStrategy diagonal() {
        return DIAGONAL;
    }

    static OUActualizationStrategy block() {
        return BLOCK;
    }

    public static OUStrategyBundle diagonalBundle() {
        return new OUStrategyBundle(OUStrategyBundle.Kind.DIAGONAL, DIAGONAL_GRADIENT);
    }

    public static OUStrategyBundle generalBundle() {
        return new OUStrategyBundle(OUStrategyBundle.Kind.GENERAL, GENERAL_GRADIENT);
    }

    public static OUStrategyBundle blockBundle() {
        return new OUStrategyBundle(OUStrategyBundle.Kind.BLOCK, BLOCK_GRADIENT);
    }

    public static OUStrategyBundle bundleFor(final MultivariateElasticModel elasticModel) {
        if (elasticModel.isDiagonal()) {
            return diagonalBundle();
        }
        if (elasticModel.hasBlockStructure() && elasticModel.hasOrthogonalActualizationBasis()) {
            return blockBundle();
        }
        return generalBundle();
    }
}
