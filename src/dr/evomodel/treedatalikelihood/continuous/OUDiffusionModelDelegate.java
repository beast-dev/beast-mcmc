/*
 * OUDiffusionModelDelegate.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 *
 */

package dr.evomodel.treedatalikelihood.continuous;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.continuous.MultivariateElasticModel;
import dr.evomodel.treedatalikelihood.continuous.cdi.ContinuousDiffusionIntegrator;
import dr.evomodel.treedatalikelihood.continuous.cdi.MultivariateIntegrator;
import dr.evomodel.treedatalikelihood.continuous.cdi.OUActualizationStrategies;
import dr.evomodel.treedatalikelihood.continuous.cdi.OUStrategyBundle;
import dr.evomodel.treedatalikelihood.continuous.cdi.SafeMultivariateDiagonalActualizedWithDriftIntegrator;
import dr.evomodel.treedatalikelihood.continuous.cdi.SafeMultivariateActualizedWithDriftIntegrator;
import dr.evomodel.treedatalikelihood.preorder.BranchSufficientStatistics;
import dr.inference.model.AbstractBlockDiagonalTwoByTwoMatrixParameter;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalLocalTransitionAdjoints;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.math.matrixAlgebra.missingData.MissingOps;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import java.util.List;

import static dr.math.matrixAlgebra.missingData.MissingOps.wrap;

/**
 * A simple OU diffusion model delegate with branch-specific drift and constant diffusion
 *
 * @author Marc A. Suchard
 * @author Paul Bastide
 */
public class OUDiffusionModelDelegate extends AbstractDriftDiffusionModelDelegate {

    private static final String USE_KERNEL_BRIDGE_PROPERTY =
            "dr.evomodel.treedatalikelihood.continuous.useTimeSeriesOUBridge";
    private static final String VALIDATE_KERNEL_BRIDGE_PROPERTY =
            "dr.evomodel.treedatalikelihood.continuous.validateTimeSeriesOUBridge";
    private static final double KERNEL_BRIDGE_TOLERANCE = 1e-8;

    // Here, branchRateModels represents optimal values

    private MultivariateElasticModel elasticModel;
    private final GaussianBranchTransitionProvider branchTransitionProvider;
    private final TimeSeriesOUCanonicalBranchGradientBridge canonicalBranchGradientBridge;

    public OUDiffusionModelDelegate(Tree tree,
                                    MultivariateDiffusionModel diffusionModel,
                                    List<BranchRateModel> branchRateModels,
                                    MultivariateElasticModel elasticModel) {
        this(tree, diffusionModel, branchRateModels, elasticModel, 0);
    }

    private OUDiffusionModelDelegate(Tree tree,
                                     MultivariateDiffusionModel diffusionModel,
                                     List<BranchRateModel> branchRateModels,
                                     MultivariateElasticModel elasticModel,
                                     int partitionNumber) {
        super(tree, diffusionModel, branchRateModels, partitionNumber);
        this.elasticModel = elasticModel;
        this.branchTransitionProvider =
                new TimeSeriesOUGaussianBranchTransitionProvider(elasticModel, diffusionModel);
        this.canonicalBranchGradientBridge =
                new TimeSeriesOUCanonicalBranchGradientBridge(
                        (TimeSeriesOUGaussianBranchTransitionProvider) this.branchTransitionProvider);
        addModel(elasticModel);
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == elasticModel) {
            fireModelChanged(model);
        } else {
            super.handleModelChangedEvent(model, object, index);
        }
    }

    @Override
    public boolean hasDrift() {
        return true;
    }

    @Override
    public boolean hasActualization() {
        return true;
    }

    @Override
    public boolean hasDiagonalActualization() {
        return elasticModel.isDiagonal();
    }

    public boolean isSymmetric() {
        return elasticModel.isSymmetric();
    }

    public double[][] getStrengthOfSelection() {
        return elasticModel.getStrengthOfSelectionMatrix();
    }

    public double[] getEigenValuesStrengthOfSelection() {
        return elasticModel.getEigenValuesStrengthOfSelection();
    }

    public double[] getEigenVectorsStrengthOfSelection() {
        return elasticModel.getEigenVectorsStrengthOfSelection();
    }

    public double[] getActualizationBasisValues() {
        return elasticModel.getActualizationBasisValues();
    }

    public double[] getActualizationBasisRotations() {
        return elasticModel.getActualizationBasisRotations();
    }

    public GaussianBranchTransitionProvider getBranchTransitionProvider() {
        return branchTransitionProvider;
    }

    public OUStrategyBundle getOUStrategyBundle() {
        return OUActualizationStrategies.bundleFor(elasticModel);
    }

    // Compatibility accessors for the legacy OU strategy helpers in cdi/.
    public Tree getTree() {
        return tree;
    }

    public MultivariateElasticModel getElasticModel() {
        return elasticModel;
    }

    public Parameter getCanonicalStationaryMeanParameter() {
        return getCanonicalConstantDriftParameter();
    }

    public int getDimTrait() {
        return dim;
    }

    public boolean hasFixedRoot(final ContinuousDataLikelihoodDelegate likelihoodDelegate) {
        return likelihoodDelegate.getRootProcessDelegate().getPseudoObservations() == Double.POSITIVE_INFINITY;
    }

    public boolean usesDiagonalActualizationBuffer(final ContinuousDiffusionIntegrator cdi) {
        return cdi instanceof SafeMultivariateDiagonalActualizedWithDriftIntegrator;
    }

    public DenseMatrix64F rootGradientVarianceWrtVariance(final NodeRef node,
                                                          final ContinuousDiffusionIntegrator cdi,
                                                          final ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                                          final DenseMatrix64F gradient) {
        return super.getGradientVarianceWrtVariance(node, cdi, likelihoodDelegate, gradient);
    }

    public double[] getBranchDriftRate(final NodeRef node) {
        return getDriftRate(node);
    }

    public boolean isExternalNode(final NodeRef node) {
        return tree.isExternal(node);
    }

    public boolean usesCanonicalVarianceGradient() {
        return elasticModel.hasBlockStructure();
    }

    public double[] getGradientVarianceForBranch(final BranchSufficientStatistics statistics,
                                                 final NodeRef node,
                                                 final ContinuousDiffusionIntegrator cdi,
                                                 final ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                                 final DenseMatrix64F gradQInv) {
        if (usesCanonicalVarianceGradient()) {
            return getCanonicalGradientVarianceForBranch(statistics, node, cdi);
        }
        return getGradientVarianceWrtVariance(node, cdi, likelihoodDelegate, gradQInv).getData();
    }

    public double[] getCanonicalGradientVarianceForBranch(final BranchSufficientStatistics statistics,
                                                          final NodeRef node,
                                                          final ContinuousDiffusionIntegrator cdi) {
        if (tree.isRoot(node)) {
            return new double[dim * dim];
        }
        final double branchLength = cdi.getBranchLength(getMatrixBufferOffsetIndex(node.getNumber()));
        final double[] optimum = getDriftRate(node);
        return canonicalBranchGradientBridge.getGradientWrtVariance(branchLength, optimum, statistics);
    }

    public double[] getCanonicalGradientDisplacementForBranch(final BranchSufficientStatistics statistics) {
        return canonicalBranchGradientBridge.getGradientWrtBranchDrift(statistics);
    }

    public double[] getCanonicalGradientDisplacementForBranch(final BranchSufficientStatistics statistics,
                                                              final NodeRef node,
                                                              final ContinuousDiffusionIntegrator cdi) {
        final double branchLength = cdi.getBranchLength(getMatrixBufferOffsetIndex(node.getNumber()));
        final double[] optimum = getDriftRate(node);
        return canonicalBranchGradientBridge.getGradientWrtBranchDrift(
                branchLength,
                optimum,
                statistics,
                node.getNumber());
    }

    public void fillCanonicalLocalAdjointsForBranch(final BranchSufficientStatistics statistics,
                                                    final CanonicalLocalTransitionAdjoints out) {
        canonicalBranchGradientBridge.fillLocalAdjoints(statistics, out);
    }

    public double[] getCanonicalGradientSelectionForBranch(final BranchSufficientStatistics statistics,
                                                           final NodeRef node,
                                                           final ContinuousDiffusionIntegrator cdi) {
        final double branchLength = cdi.getBranchLength(getMatrixBufferOffsetIndex(node.getNumber()));
        final double[] optimum = getDriftRate(node);
        return canonicalBranchGradientBridge.getGradientWrtSelection(branchLength, optimum, statistics);
    }

    public double[] getCanonicalNativeOrthogonalBlockGradientSelectionForBranch(final BranchSufficientStatistics statistics,
                                                                                final NodeRef node,
                                                                                final ContinuousDiffusionIntegrator cdi,
                                                                                final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter,
                                                                                final Parameter requestedParameter) {
        final double branchLength = cdi.getBranchLength(getMatrixBufferOffsetIndex(node.getNumber()));
        final double[] optimum = getDriftRate(node);
        return canonicalBranchGradientBridge.getGradientWrtSelection(
                branchLength,
                optimum,
                statistics,
                tree.isExternal(node),
                blockParameter,
                requestedParameter);
    }

    /**
     * Accumulates the global pruning-remainder selection gradient for one branch.
     */
    public void accumulateCanonicalGlobalRemainderSelectionGradientForBranch(
            final NodeRef node,
            final NodeRef parentNode,
            final int nodeBuffer,
            final int parentBuffer,
            final ContinuousDiffusionIntegrator cdi,
            final AbstractBlockDiagonalTwoByTwoMatrixParameter nativeBlockParameter,
            final Parameter requestedParameter,
            final double[] gradOut) {

        final int matrixIndex = getMatrixBufferOffsetIndex(node.getNumber());
        final double branchLength = cdi.getBranchLength(matrixIndex);
        final int dimTrait = canonicalBranchGradientBridge.getDimension();
        final double[] barVdi = new double[dimTrait * dimTrait];

        cdi.computeGlobalRemainderAdjointWrtBranchVariance(
                nodeBuffer, matrixIndex, parentBuffer, barVdi);

        canonicalBranchGradientBridge.accumulateGlobalRemainderSelectionGradientForBranch(
                branchLength, barVdi, nativeBlockParameter, requestedParameter, gradOut);
    }

    @Override
    public void setDiffusionModels(ContinuousDiffusionIntegrator cdi, boolean flip) {
        super.setDiffusionModels(cdi, flip);

        cdi.setDiffusionStationaryVariance(getEigenBufferOffsetIndex(0),
                getActualizationBasisValues(), getActualizationBasisRotations());
    }

    @Override
    public void updateDiffusionMatrices(ContinuousDiffusionIntegrator cdi, int[] branchIndices, double[] edgeLengths,
                                        int updateCount, boolean flip) {

        int[] probabilityIndices = new int[updateCount];
        final double[] driftRates = getDriftRates(branchIndices, updateCount);

        for (int i = 0; i < updateCount; i++) {
            if (flip) {
                flipMatrixBufferOffset(branchIndices[i]);
            }
            probabilityIndices[i] = getMatrixBufferOffsetIndex(branchIndices[i]);
        }

        cdi.updateOrnsteinUhlenbeckDiffusionMatrices(
                getEigenBufferOffsetIndex(0),
                probabilityIndices,
                edgeLengths,
                driftRates,
                getActualizationBasisValues(),
                getActualizationBasisRotations(),
                updateCount);

        if (Boolean.getBoolean(USE_KERNEL_BRIDGE_PROPERTY)) {
            applyKernelBridgeForwardPath(cdi, probabilityIndices, edgeLengths, driftRates, updateCount);
        }

        if (Boolean.getBoolean(VALIDATE_KERNEL_BRIDGE_PROPERTY)) {
            validateKernelBridge(cdi, branchIndices, probabilityIndices, edgeLengths, updateCount);
        }
    }

    private void applyKernelBridgeForwardPath(final ContinuousDiffusionIntegrator cdi,
                                              final int[] probabilityIndices,
                                              final double[] edgeLengths,
                                              final double[] optimalRates,
                                              final int updateCount) {
        if (hasDiagonalActualization()) {
            if (!(cdi instanceof SafeMultivariateDiagonalActualizedWithDriftIntegrator)) {
                throw new IllegalStateException("Kernel bridge replacement requires a diagonal actualized OU integrator");
            }
            applyKernelBridgeForwardPathDiagonal(
                    (SafeMultivariateDiagonalActualizedWithDriftIntegrator) cdi,
                    probabilityIndices,
                    edgeLengths,
                    optimalRates,
                    updateCount);
            return;
        }
        if (!(cdi instanceof SafeMultivariateActualizedWithDriftIntegrator)) {
            throw new IllegalStateException("Kernel bridge replacement requires a full actualized OU integrator");
        }
        applyKernelBridgeForwardPathFull(
                (SafeMultivariateActualizedWithDriftIntegrator) cdi,
                probabilityIndices,
                edgeLengths,
                optimalRates,
                updateCount);
    }

    private void applyKernelBridgeForwardPathDiagonal(
            final SafeMultivariateDiagonalActualizedWithDriftIntegrator cdi,
            final int[] probabilityIndices,
            final double[] edgeLengths,
            final double[] optimalRates,
            final int updateCount) {
        final double[][] actualizationMatrix = new double[dim][dim];
        final double[][] varianceMatrix = new double[dim][dim];
        final double[] diagonalActualization = new double[dim];
        final double[] displacement = new double[dim];
        final double[] varianceFlat = new double[dim * dim];

        for (int up = 0; up < updateCount; ++up) {
            final int probabilityIndex = probabilityIndices[up];
            final int optimalOffset = up * dim;

            branchTransitionProvider.fillBranchTransitionMatrix(edgeLengths[up], actualizationMatrix);
            branchTransitionProvider.fillBranchTransitionCovariance(edgeLengths[up], varianceMatrix);

            for (int i = 0; i < dim; ++i) {
                final double actualization = actualizationMatrix[i][i];
                diagonalActualization[i] = actualization;
                displacement[i] = (1.0 - actualization) * optimalRates[optimalOffset + i];
                for (int j = 0; j < dim; ++j) {
                    varianceFlat[i * dim + j] = varianceMatrix[i][j];
                }
            }

            cdi.overrideOrnsteinUhlenbeckBranchMatrices(
                    probabilityIndex,
                    diagonalActualization,
                    varianceFlat,
                    displacement);
        }
    }

    private void applyKernelBridgeForwardPathFull(
            final SafeMultivariateActualizedWithDriftIntegrator cdi,
            final int[] probabilityIndices,
            final double[] edgeLengths,
            final double[] optimalRates,
            final int updateCount) {
        final double[][] actualizationMatrix = new double[dim][dim];
        final double[][] varianceMatrix = new double[dim][dim];
        final double[] actualizationFlat = new double[dim * dim];
        final double[] displacement = new double[dim];
        final double[] varianceFlat = new double[dim * dim];

        for (int up = 0; up < updateCount; ++up) {
            final int probabilityIndex = probabilityIndices[up];
            final int optimalOffset = up * dim;

            branchTransitionProvider.fillBranchTransitionMatrix(edgeLengths[up], actualizationMatrix);
            branchTransitionProvider.fillBranchTransitionCovariance(edgeLengths[up], varianceMatrix);

            for (int i = 0; i < dim; ++i) {
                double propagatedOptimum = 0.0;
                for (int j = 0; j < dim; ++j) {
                    final double actualization = actualizationMatrix[i][j];
                    actualizationFlat[i * dim + j] = actualization;
                    varianceFlat[i * dim + j] = varianceMatrix[i][j];
                    propagatedOptimum += actualization * optimalRates[optimalOffset + j];
                }
                displacement[i] = optimalRates[optimalOffset + i] - propagatedOptimum;
            }

            cdi.overrideOrnsteinUhlenbeckBranchMatrices(
                    probabilityIndex,
                    actualizationFlat,
                    varianceFlat,
                    displacement);
        }
    }

    private void validateKernelBridge(final ContinuousDiffusionIntegrator cdi,
                                      final int[] branchIndices,
                                      final int[] probabilityIndices,
                                      final double[] edgeLengths,
                                      final int updateCount) {
        final double[][] expectedActualization = new double[dim][dim];
        final double[][] expectedVariance = new double[dim][dim];
        final double[] observedActualizationRaw = hasDiagonalActualization()
                ? new double[dim]
                : new double[dim * dim];
        final double[] observedVarianceRaw = new double[dim * dim];

        for (int i = 0; i < updateCount; ++i) {
            final int branchIndex = branchIndices[i];
            final int probabilityIndex = probabilityIndices[i];
            final double branchLength = edgeLengths[i];

            branchTransitionProvider.fillBranchTransitionMatrix(branchLength, expectedActualization);
            branchTransitionProvider.fillBranchTransitionCovariance(branchLength, expectedVariance);

            cdi.getBranchActualization(probabilityIndex, observedActualizationRaw);
            cdi.getBranchVariance(probabilityIndex, getEigenBufferOffsetIndex(0), observedVarianceRaw);

            assertActualizationMatches(branchIndex, expectedActualization, observedActualizationRaw);
            assertMatrixMatches("variance", branchIndex, expectedVariance, observedVarianceRaw);
        }
    }

    private void assertActualizationMatches(final int branchIndex,
                                            final double[][] expected,
                                            final double[] observedRaw) {
        if (hasDiagonalActualization()) {
            for (int i = 0; i < dim; ++i) {
                assertClose("actualization", branchIndex, i, i, expected[i][i], observedRaw[i]);
            }
            return;
        }
        assertMatrixMatches("actualization", branchIndex, expected, observedRaw);
    }

    private void assertMatrixMatches(final String label,
                                     final int branchIndex,
                                     final double[][] expected,
                                     final double[] observedRaw) {
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                assertClose(label, branchIndex, i, j, expected[i][j], observedRaw[i * dim + j]);
            }
        }
    }

    private static void assertClose(final String label,
                                    final int branchIndex,
                                    final int row,
                                    final int col,
                                    final double expected,
                                    final double observed) {
        final double scale = Math.max(1.0, Math.max(Math.abs(expected), Math.abs(observed)));
        if (Math.abs(expected - observed) > KERNEL_BRIDGE_TOLERANCE * scale) {
            throw new IllegalStateException(
                    "Kernel bridge mismatch for " + label +
                            " on branch " + branchIndex +
                            " at [" + row + "," + col + "]: expected=" + expected +
                            ", observed=" + observed);
        }
    }

    ///////////////////////////////////////////////////////////////////////////
    /// Gradient Functions
    ///////////////////////////////////////////////////////////////////////////

    @Override
    public DenseMatrix64F getGradientVarianceWrtVariance(NodeRef node,
                                                         ContinuousDiffusionIntegrator cdi,
                                                         ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                                         DenseMatrix64F gradient) {
        if (tree.isRoot(node)) {
            return super.getGradientVarianceWrtVariance(node, cdi, likelihoodDelegate, gradient);
        } else {
            DenseMatrix64F result = gradient.copy();
            if (hasDiagonalActualization()) {
                actualizeGradientDiagonal(cdi, node.getNumber(), result);
            } else {
                actualizeGradient(cdi, node.getNumber(), result);
            }
            return result;
        }
    }

    private void actualizeGradient(ContinuousDiffusionIntegrator cdi, int nodeIndex, DenseMatrix64F gradient) {
        double[] attenuationRotation = elasticModel.getEigenVectorsStrengthOfSelection();
        DenseMatrix64F P = wrap(attenuationRotation, 0, dim, dim);

        SafeMultivariateActualizedWithDriftIntegrator.transformMatrix(gradient, P, elasticModel.isSymmetric());

        actualizeGradientDiagonal(cdi, nodeIndex, gradient);

        SafeMultivariateActualizedWithDriftIntegrator.transformMatrixBack(gradient, P);

    }

    private void actualizeGradientDiagonal(ContinuousDiffusionIntegrator cdi, int nodeIndex, DenseMatrix64F gradient) {
        double[] attenuation = elasticModel.getEigenValuesStrengthOfSelection();
        double edgeLength = cdi.getBranchLength(getMatrixBufferOffsetIndex(nodeIndex));
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                gradient.unsafe_set(i, j,
                        factorFunction(attenuation[i] + attenuation[j], edgeLength) * gradient.unsafe_get(i, j));
            }
        }
    }

    private static double factorFunction(double x, double l) {
        if (x == 0) return l;
        return -Math.expm1(-x * l) / x;
    }

//    private void actualizeGradientOld(ContinuousDiffusionIntegrator cdi, int nodeIndex, DenseMatrix64F gradient) {
//        double[] actualization = new double[dim * dim];
//        cdi.getBranchActualization(getMatrixBufferOffsetIndex(nodeIndex), actualization);
//        double[] attenuation = elasticModel.getStrengthOfSelectionMatrixAsVector();
//        // Actualization
//        DenseMatrix64F A = DenseMatrix64F.wrap(dim, dim, actualization);
//        DenseMatrix64F temp = new DenseMatrix64F(dim, dim);
//        CommonOps.multTransB(gradient, A, temp);
//        CommonOps.multAdd(-1.0, A, temp, gradient);
//        // Derivation of stationary variance wrt variance
//        DenseMatrix64F factor = DenseMatrix64F.wrap(dim * dim, dim * dim,
//                KroneckerOperation.sum(attenuation, dim, attenuation, dim));
//        CommonOps.invert(factor);
//        DenseMatrix64F temp1 = DenseMatrix64F.wrap(dim * dim, 1, gradient.getData());
//        DenseMatrix64F temp2 = new DenseMatrix64F(dim * dim, 1);
//        CommonOps.mult(factor, temp1, temp2);
//        gradient.setData(temp2.getData());
//    }

    DenseMatrix64F getGradientVarianceWrtAttenuation(NodeRef node,
                                                     ContinuousDiffusionIntegrator cdi,
                                                     BranchSufficientStatistics statistics,
                                                     DenseMatrix64F gradient) {
        assert !tree.isRoot(node) : "Gradient wrt actualization is not available for the root.";
        if (hasDiagonalActualization()) {
            return getGradientVarianceWrtAttenuationDiagonal(cdi, statistics, node.getNumber(), gradient);
        } else {
            throw new RuntimeException("not yet implemented");
        }
    }

    private DenseMatrix64F getGradientVarianceWrtAttenuationDiagonal(ContinuousDiffusionIntegrator cdi,
                                                                     BranchSufficientStatistics statistics,
                                                                     int nodeIndex, DenseMatrix64F gradient) {
        // wrt to q_i actualization
        DenseMatrix64F gradActualization = getGradientVarianceWrtActualizationDiagonal(cdi, statistics, nodeIndex, gradient);

        // wrt to Gamma stationary variance
        DenseMatrix64F gradStationary = getGradientBranchVarianceWrtAttenuationDiagonal(cdi, nodeIndex, gradient);

        CommonOps.addEquals(gradActualization, gradStationary);
        return gradActualization;

    }

    private DenseMatrix64F getGradientVarianceWrtActualizationDiagonal(ContinuousDiffusionIntegrator cdi, BranchSufficientStatistics statistics,
                                                                       int nodeIndex, DenseMatrix64F gradient) {
        // q_i
//        double[] qi = new double[dim];
//        cdi.getBranchActualization(getMatrixBufferOffsetIndex(nodeIndex), qi);
//        DenseMatrix64F qi = wrapDiagonal(actualization, 0, dim);
        // q_i^-1
//        DenseMatrix64F qiInv = wrapDiagonalInverse(actualization, 0, dim);
        // Q_i^-
        DenseMatrix64F Wi = statistics.getAbove().getRawVarianceCopy();
        // Gamma
//        DenseMatrix64F Gamma = wrap(
//                ((SafeMultivariateDiagonalActualizedWithDriftIntegrator) cdi).getStationaryVariance(getEigenBufferOffsetIndex(0)),
//                0, dim, dim);
        // Branch variance
        double[] branchVariance = new double[dim * dim];
        cdi.getBranchVariance(getMatrixBufferOffsetIndex(nodeIndex), getEigenBufferOffsetIndex(0), branchVariance);
        DenseMatrix64F Sigma_i = wrap(branchVariance, 0, dim, dim);

        // factor
        DenseMatrix64F res = new DenseMatrix64F(dim, dim);
        CommonOps.addEquals(Wi, -1, Sigma_i);
//        MissingOps.diagDiv(qi, Wi);
//        CommonOps.multTransA(qiInv, tmp, res);
        CommonOps.multTransB(Wi, gradient, res);

        // temp + temp^T
//        MissingOps.addTransEquals(res); No need for diagonal case
        CommonOps.scale(2.0, res);

        // Get diagonal
        DenseMatrix64F resDiag = new DenseMatrix64F(dim, 1);
        CommonOps.extractDiag(res, resDiag);

        // Wrt attenuation
        double ti = cdi.getBranchLength(getMatrixBufferOffsetIndex(nodeIndex));
        chainRuleActualizationWrtAttenuationDiagonal(ti, resDiag);

        return resDiag;

    }

    private void chainRuleActualizationWrtAttenuationDiagonal(double ti, DenseMatrix64F grad) {
//        MissingOps.diagMult(actualization, grad);
        CommonOps.scale(-ti, grad);
    }

    private DenseMatrix64F getGradientBranchVarianceWrtAttenuationDiagonal(ContinuousDiffusionIntegrator cdi,
                                                                           int nodeIndex, DenseMatrix64F gradient) {

        double[] attenuation = elasticModel.getEigenValuesStrengthOfSelection();
        DenseMatrix64F variance = wrap(
                ((MultivariateIntegrator) cdi).getVariance(getEigenBufferOffsetIndex(0)),
                0, dim, dim);

        double ti = cdi.getBranchLength(getMatrixBufferOffsetIndex(nodeIndex));

        DenseMatrix64F res = new DenseMatrix64F(dim, 1);

        CommonOps.elementMult(variance, gradient);

        for (int k = 0; k < dim; k++) {
            double sum = 0.0;
            for (int l = 0; l < dim; l++) {
                sum -= variance.unsafe_get(k, l) * computeAttenuationFactorActualized(attenuation[k] + attenuation[l], ti);
            }
            res.unsafe_set(k, 0, sum);
        }

        return res;
    }

    private double computeAttenuationFactorActualized(double lambda, double ti) {
        if (lambda == 0) return ti * ti;
        double em1 = Math.expm1(-lambda * ti);
        return 2.0 * (em1 * em1 - (em1 + lambda * ti) * Math.exp(-lambda * ti)) / lambda / lambda;
//        return 2.0 * (1 - (1 + lambda * ti) * Math.exp( - lambda * ti)) / (lambda * lambda);
    }

    DenseMatrix64F getGradientDisplacementWrtAttenuation(NodeRef node,
                                                         ContinuousDiffusionIntegrator cdi,
                                                         BranchSufficientStatistics statistics,
                                                         DenseMatrix64F gradient) {
        assert !tree.isRoot(node) : "Gradient wrt actualization is not available for the root.";
        if (hasDiagonalActualization()) {
            return getGradientDisplacementWrtAttenuationDiagonal(cdi, statistics, node, gradient);
        } else {
            throw new RuntimeException("not yet implemented");
        }
    }

    private DenseMatrix64F getGradientDisplacementWrtAttenuationDiagonal(ContinuousDiffusionIntegrator cdi, BranchSufficientStatistics statistics,
                                                                         NodeRef node, DenseMatrix64F gradient) {
        int nodeIndex = node.getNumber();
        // q_i
//        double[] qi = new double[dim];
//        cdi.getBranchActualization(getMatrixBufferOffsetIndex(nodeIndex), qi);
//        DenseMatrix64F qi = wrapDiagonal(actualization, 0, dim);
        // q_i^-1
//        DenseMatrix64F qiInv = wrapDiagonalInverse(actualization, 0, dim);
        // ni
        DenseMatrix64F ni = statistics.getAbove().getRawMean();
        // beta_i
        DenseMatrix64F betai = wrap(getDriftRate(node), 0, dim, 1);

        // factor
//        DenseMatrix64F tmp = new DenseMatrix64F(dim, 1);
        DenseMatrix64F resFull = new DenseMatrix64F(dim, dim);
        DenseMatrix64F resDiag = new DenseMatrix64F(dim, 1);
        CommonOps.add(ni, -1, betai, resDiag);
//        MissingOps.diagDiv(qi, resDiag);
        CommonOps.multTransB(gradient, resDiag, resFull);

        // Extract diag
        CommonOps.extractDiag(resFull, resDiag);

        // Wrt attenuation
        double ti = cdi.getBranchLength(getMatrixBufferOffsetIndex(nodeIndex));
        chainRuleActualizationWrtAttenuationDiagonal(ti, resDiag);

        return resDiag;

    }

    @Override
    DenseMatrix64F getGradientDisplacementWrtDrift(NodeRef node,
                                                   ContinuousDiffusionIntegrator cdi,
                                                   ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                                   DenseMatrix64F gradient) {
        DenseMatrix64F result = gradient.copy();
        if (hasDiagonalActualization()) {
            actualizeDisplacementGradientDiagonal(cdi, node.getNumber(), result);
        } else {
            actualizeDisplacementGradient(cdi, node.getNumber(), result);
        }
        return result;
    }

    private void actualizeDisplacementGradientDiagonal(ContinuousDiffusionIntegrator cdi,
                                                       int nodeIndex, DenseMatrix64F gradient) {
        // q_i
        double[] qi = new double[dim];
        cdi.getBranch1mActualization(getMatrixBufferOffsetIndex(nodeIndex), qi);
        MissingOps.diagMult(qi, gradient);
    }

    private void actualizeDisplacementGradient(ContinuousDiffusionIntegrator cdi,
                                               int nodeIndex, DenseMatrix64F gradient) {
        // q_i
        double[] qi = new double[dim * dim];
        cdi.getBranch1mActualization(getMatrixBufferOffsetIndex(nodeIndex), qi);
        DenseMatrix64F Actu = wrap(qi, 0, dim, dim);
        CommonOps.scale(-1.0, Actu);
//        for (int i = 0; i < dim; i++) {
//            Actu.unsafe_set(i, i, Actu.unsafe_get(i, i) - 1.0);
//        }
        DenseMatrix64F tmp = new DenseMatrix64F(dim, 1);
        CommonOps.mult(Actu, gradient, tmp);
        CommonOps.scale(-1.0, tmp, gradient);
    }

    @Override
    public double[] getGradientDisplacementWrtRoot(NodeRef node,
                                                   ContinuousDiffusionIntegrator cdi,
                                                   ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                                   DenseMatrix64F gradient) {
        boolean fixedRoot = likelihoodDelegate.getRootProcessDelegate().getPseudoObservations() == Double.POSITIVE_INFINITY;
        if (fixedRoot && tree.isRoot(tree.getParent(node))) {
            return actualizeRootGradient(cdi, node.getNumber(), gradient);
        }
        if (!fixedRoot && tree.isRoot(node)) {
            return gradient.getData();
        }
        return new double[gradient.getNumRows()];
    }

    private double[] actualizeRootGradient(ContinuousDiffusionIntegrator cdi,
                                           int nodeIndex,
                                           DenseMatrix64F gradient) {
        if (hasDiagonalActualization()) {
            return actualizeRootGradientDiagonal(cdi, nodeIndex, gradient);
        } else {
            return actualizeRootGradientFull(cdi, nodeIndex, gradient);
        }
    }

    private double[] actualizeRootGradientDiagonal(ContinuousDiffusionIntegrator cdi,
                                                   int nodeIndex, DenseMatrix64F gradient) {
        double[] qi = new double[dim];
        cdi.getBranchActualization(getMatrixBufferOffsetIndex(nodeIndex), qi);
        DenseMatrix64F tmp = new DenseMatrix64F(dim, 1);
        MissingOps.diagMult(qi, gradient, tmp);
        return tmp.getData();
    }

    private double[] actualizeRootGradientFull(ContinuousDiffusionIntegrator cdi,
                                               int nodeIndex, DenseMatrix64F gradient) {
        // q_i
        double[] qi = new double[dim * dim];
        cdi.getBranchActualization(getMatrixBufferOffsetIndex(nodeIndex), qi);
        DenseMatrix64F Actu = wrap(qi, 0, dim, dim);
        DenseMatrix64F tmp = new DenseMatrix64F(dim, 1);
        // Root-mean backprop is a vector-Jacobian product; use A^T * g for full (non-symmetric) actualization.
        CommonOps.multTransA(Actu, gradient, tmp);
        return tmp.getData();
    }


    ///////////////////////////////////////////////////////////////////////////
    /// Report Functions
    ///////////////////////////////////////////////////////////////////////////

    @Override
    public double[][] getJointVariance(final double priorSampleSize,
                                       final double[][] treeVariance, final double[][] treeSharedLengths,
                                       final double[][] traitVariance) {
        if (hasDiagonalActualization()) {
            return getJointVarianceDiagonal(priorSampleSize, treeVariance, treeSharedLengths, traitVariance);
        }
        return getJointVarianceFull(priorSampleSize, treeVariance, treeSharedLengths, traitVariance);
    }

    private double[][] getJointVarianceFull(final double priorSampleSize,
                                            final double[][] treeVariance, final double[][] treeSharedLengths,
                                            final double[][] traitVariance) {

        double[] eigVals = this.getEigenValuesStrengthOfSelection();
        DenseMatrix64F V = wrap(this.getEigenVectorsStrengthOfSelection(), 0, dim, dim);
        DenseMatrix64F Vinv = new DenseMatrix64F(dim, dim);
        CommonOps.invert(V, Vinv);

        DenseMatrix64F transTraitVariance = new DenseMatrix64F(traitVariance);

        DenseMatrix64F tmp = new DenseMatrix64F(dim, dim);
        CommonOps.mult(Vinv, transTraitVariance, tmp);
        CommonOps.multTransB(tmp, Vinv, transTraitVariance);

        // inverse of eigenvalues
        double[][] invEigVals = new double[dim][dim];
        for (int p = 0; p < dim; ++p) {
            for (int q = 0; q < dim; ++q) {
                invEigVals[p][q] = 1 / (eigVals[p] + eigVals[q]);
            }
        }

        // Computation of matrix
        int ntaxa = tree.getExternalNodeCount();
        double ti;
        double tj;
        double tij;
        double ep;
        double eq;
        DenseMatrix64F varTemp = new DenseMatrix64F(dim, dim);
        double[][] jointVariance = new double[dim * ntaxa][dim * ntaxa];
        for (int i = 0; i < ntaxa; ++i) {
            for (int j = 0; j < ntaxa; ++j) {
                ti = treeSharedLengths[i][i];
                tj = treeSharedLengths[j][j];
                tij = treeSharedLengths[i][j];
                for (int p = 0; p < dim; ++p) {
                    for (int q = 0; q < dim; ++q) {
                        ep = eigVals[p];
                        eq = eigVals[q];
                        varTemp.set(p, q, Math.exp(-ep * ti) * Math.exp(-eq * tj) * (invEigVals[p][q] * (Math.exp((ep + eq) * tij) - 1) + 1 / priorSampleSize) * transTraitVariance.get(p, q));
                    }
                }
                CommonOps.mult(V, varTemp, tmp);
                CommonOps.multTransB(tmp, V, varTemp);
                for (int p = 0; p < dim; ++p) {
                    for (int q = 0; q < dim; ++q) {
                        jointVariance[i * dim + p][j * dim + q] = varTemp.get(p, q);
                    }
                }
            }
        }
        return jointVariance;
    }

    private double[][] getJointVarianceDiagonal(final double priorSampleSize,
                                                final double[][] treeVariance, final double[][] treeSharedLengths,
                                                final double[][] traitVariance) {

        // Eigen of strength of selection matrix
        double[] eigVals = this.getEigenValuesStrengthOfSelection();
        int ntaxa = tree.getExternalNodeCount();
        double ti;
        double tj;
        double tij;
        double ep;
        double eq;
        double var;
        DenseMatrix64F varTemp = new DenseMatrix64F(dim, dim);
        double[][] jointVariance = new double[dim * ntaxa][dim * ntaxa];
        for (int i = 0; i < ntaxa; ++i) {
            for (int j = 0; j < ntaxa; ++j) {
                ti = treeSharedLengths[i][i];
                tj = treeSharedLengths[j][j];
                tij = treeSharedLengths[i][j];
                for (int p = 0; p < dim; ++p) {
                    for (int q = 0; q < dim; ++q) {
                        ep = eigVals[p];
                        eq = eigVals[q];
                        if (ep + eq == 0.0) {
                            var = (tij + 1 / priorSampleSize) * traitVariance[p][q];
                        } else {
                            var = Math.exp(-ep * ti) * Math.exp(-eq * tj) * (Math.expm1((ep + eq) * tij) / (ep + eq) + 1 / priorSampleSize) * traitVariance[p][q];
                        }
                        varTemp.set(p, q, var);
                    }
                }
                for (int p = 0; p < dim; ++p) {
                    for (int q = 0; q < dim; ++q) {
                        jointVariance[i * dim + p][j * dim + q] = varTemp.get(p, q);
                    }
                }
            }
        }
        return jointVariance;
    }

    ///////////////////////////////////////////////////////////////////////////
    /// Model Heritability Functions
    ///////////////////////////////////////////////////////////////////////////

    @Override
    public void getMeanTipVariances(final double priorSampleSize,
                                    final double[] treeLengths,
                                    final DenseMatrix64F traitVariance,
                                    final DenseMatrix64F varSum) {
        if (hasDiagonalActualization()) {
            getMeanTipVariancesDiagonal(priorSampleSize, treeLengths, traitVariance, varSum);
        }
        getMeanTipVariancesFull(priorSampleSize, treeLengths, traitVariance, varSum);
    }

    private void getMeanTipVariancesFull(final double priorSampleSize,
                                         final double[] treeLengths,
                                         final DenseMatrix64F traitVariance,
                                         final DenseMatrix64F varSum) {

        DenseMatrix64F V = wrap(this.getEigenVectorsStrengthOfSelection(), 0, dim, dim);
        DenseMatrix64F Vinv = new DenseMatrix64F(dim, dim);
        CommonOps.invert(V, Vinv);

        DenseMatrix64F transTraitVariance = new DenseMatrix64F(traitVariance);

        DenseMatrix64F tmp = new DenseMatrix64F(dim, dim);
        CommonOps.mult(Vinv, transTraitVariance, tmp);
        CommonOps.multTransB(tmp, Vinv, transTraitVariance);

        // Diagonal Computations
        getMeanTipVariancesDiagonal(priorSampleSize, treeLengths, transTraitVariance, varSum);

        // Back to original space
        CommonOps.mult(V, varSum, tmp);
        CommonOps.multTransB(tmp, V, varSum);
    }

    private void getMeanTipVariancesDiagonal(final double priorSampleSize,
                                             final double[] treeLengths,
                                             final DenseMatrix64F traitVariance,
                                             final DenseMatrix64F varSum) {

        // Eigen of strength of selection matrix
        double[] eigVals = this.getEigenValuesStrengthOfSelection();
        int ntaxa = tree.getExternalNodeCount();
        double ti;
        double ep;
        double eq;
        double var;

        for (int i = 0; i < ntaxa; ++i) {
            ti = treeLengths[i];
            for (int p = 0; p < dim; ++p) {
                ep = eigVals[p];
                for (int q = 0; q < dim; ++q) {
                    eq = eigVals[q];
                    if (ep + eq == 0.0) {
                        var = (ti + 1 / priorSampleSize) * traitVariance.get(p, q);
                    } else {
                        var = Math.exp(-(ep + eq) * ti) * (Math.expm1((ep + eq) * ti) / (ep + eq) + 1 / priorSampleSize) * traitVariance.get(p, q);
                    }
                    varSum.set(p, q, varSum.get(p, q) + var);
                }
            }
        }
        CommonOps.scale(1.0 / treeLengths.length, varSum);
    }
}
