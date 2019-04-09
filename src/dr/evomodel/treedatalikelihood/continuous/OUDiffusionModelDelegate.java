/*
 * AbstractOUDiffusionModelDelegate.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
 */

package dr.evomodel.treedatalikelihood.continuous;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.continuous.MultivariateElasticModel;
import dr.evomodel.treedatalikelihood.continuous.cdi.ContinuousDiffusionIntegrator;
import dr.evomodel.treedatalikelihood.continuous.cdi.MultivariateIntegrator;
import dr.evomodel.treedatalikelihood.continuous.cdi.SafeMultivariateActualizedWithDriftIntegrator;
import dr.evomodel.treedatalikelihood.continuous.cdi.SafeMultivariateDiagonalActualizedWithDriftIntegrator;
import dr.evomodel.treedatalikelihood.preorder.BranchSufficientStatistics;
import dr.inference.model.Model;
import dr.math.KroneckerOperation;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import java.util.List;

import static dr.math.matrixAlgebra.missingData.MissingOps.wrap;
import static dr.math.matrixAlgebra.missingData.MissingOps.wrapDiagonal;

/**
 * A simple OU diffusion model delegate with branch-specific drift and constant diffusion
 *
 * @author Marc A. Suchard
 * @author Paul Bastide
 * @version $Id$
 */
public class OUDiffusionModelDelegate extends AbstractDriftDiffusionModelDelegate {

    // Here, branchRateModels represents optimal values

    private MultivariateElasticModel elasticModel;

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

    @Override
    public void setDiffusionModels(ContinuousDiffusionIntegrator cdi, boolean flip) {
        super.setDiffusionModels(cdi, flip);

        cdi.setDiffusionStationaryVariance(getEigenBufferOffsetIndex(0),
                getEigenValuesStrengthOfSelection(), getEigenVectorsStrengthOfSelection());
    }

    @Override
    public void updateDiffusionMatrices(ContinuousDiffusionIntegrator cdi, int[] branchIndices, double[] edgeLengths,
                                        int updateCount, boolean flip) {

        int[] probabilityIndices = new int[updateCount];

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
                getDriftRates(branchIndices, updateCount),
                getEigenValuesStrengthOfSelection(),
                getEigenVectorsStrengthOfSelection(),
                updateCount);
    }

    ///////////////////////////////////////////////////////////////////////////
    /// Gradient Functions
    ///////////////////////////////////////////////////////////////////////////

    @Override
    public void getGradientVarianceWrtVariance(NodeRef node,
                                               ContinuousDiffusionIntegrator cdi,
                                               ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                               DenseMatrix64F gradient) {
        if (tree.isRoot(node)) {
            super.getGradientVarianceWrtVariance(node, cdi, likelihoodDelegate, gradient);
        } else {
            if (hasDiagonalActualization()) {
                actualizeGradientDiagonal(cdi, node.getNumber(), gradient);
            } else {
                actualizeGradient(cdi, node.getNumber(), gradient);
            }
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
        return (1 - Math.exp(-x * l)) / x;
    }

    private void actualizeGradientOld(ContinuousDiffusionIntegrator cdi, int nodeIndex, DenseMatrix64F gradient) {
        double[] actualization = new double[dim * dim];
        cdi.getBranchActualization(getMatrixBufferOffsetIndex(nodeIndex), actualization);
        double[] attenuation = elasticModel.getStrengthOfSelectionMatrixAsVector();
        // Actualization
        DenseMatrix64F A = DenseMatrix64F.wrap(dim, dim, actualization);
        DenseMatrix64F temp = new DenseMatrix64F(dim, dim);
        CommonOps.multTransB(gradient, A, temp);
        CommonOps.multAdd(-1.0, A, temp, gradient);
        // Derivation of stationary variance wrt variance
        DenseMatrix64F factor = DenseMatrix64F.wrap(dim * dim, dim * dim,
                KroneckerOperation.sum(attenuation, dim, attenuation, dim));
        CommonOps.invert(factor);
        DenseMatrix64F temp1 = DenseMatrix64F.wrap(dim * dim, 1, gradient.getData());
        DenseMatrix64F temp2 = new DenseMatrix64F(dim * dim, 1);
        CommonOps.mult(factor, temp1, temp2);
        gradient.setData(temp2.getData());
    }

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
        DenseMatrix64F gradStationary = getGradientVarianceWrtStationaryDiagonal(cdi, nodeIndex, gradient);

        DenseMatrix64F res = new DenseMatrix64F(dim, 1);
        CommonOps.add(gradActualization, gradStationary, res);
        return (res);

    }

    private DenseMatrix64F getGradientVarianceWrtActualizationDiagonal(ContinuousDiffusionIntegrator cdi, BranchSufficientStatistics statistics,
                                                                       int nodeIndex, DenseMatrix64F gradient) {
        // q_i
        double[] actualization = new double[dim];
        cdi.getBranchActualization(getMatrixBufferOffsetIndex(nodeIndex), actualization);
        DenseMatrix64F qi = wrapDiagonal(actualization, 0, dim);
        // q_i^-1
        double[] actualizationInv = new double[dim];
        for (int i = 0; i < dim; i++) {
            actualizationInv[i] = 1.0 / actualization[i];
        }
        DenseMatrix64F qiInv = wrapDiagonal(actualizationInv, 0, dim);
        // Q_i^-
        DenseMatrix64F Wi = statistics.getAbove().getRawVariance();
        // Gamma
        DenseMatrix64F Gamma = wrap(
                ((SafeMultivariateDiagonalActualizedWithDriftIntegrator) cdi).getStationaryVariance(getEigenBufferOffsetIndex(0)),
                0, dim, dim);

        // factor
        DenseMatrix64F tmp = new DenseMatrix64F(dim, dim);
        DenseMatrix64F res = new DenseMatrix64F(dim, dim);
        CommonOps.add(Wi, -1, Gamma, tmp);
        CommonOps.mult(qiInv, tmp, res);
        CommonOps.mult(res, gradient, tmp);

        // temp + temp^T
        for (int i = 0; i < gradient.getNumCols(); i++) {
            for (int j = 0; j < gradient.getNumCols(); j++) {
                res.unsafe_set(i, j, tmp.unsafe_get(i, j) + tmp.unsafe_get(j, i));
            }
        }

        // Get diagonal
        DenseMatrix64F resDiag = new DenseMatrix64F(dim, 1);
        for (int i = 0; i < dim; ++i) {
            resDiag.unsafe_set(i, 0, res.unsafe_get(i, i));
        }

        // Wrt attenuation
        double ti = cdi.getBranchLength(getMatrixBufferOffsetIndex(nodeIndex));
        chainRuleActualizationWrtAttenuationDiagonal(qi, ti, resDiag);

        return resDiag;

    }

    private void chainRuleActualizationWrtAttenuationDiagonal(DenseMatrix64F qi, double ti,
                                                              DenseMatrix64F grad) {
        DenseMatrix64F tmp = new DenseMatrix64F(dim, 1);
        CommonOps.mult(qi, grad, tmp);
        CommonOps.scale(-ti, tmp, grad);
    }

    private DenseMatrix64F getGradientVarianceWrtStationaryDiagonal(ContinuousDiffusionIntegrator cdi,
                                                                    int nodeIndex, DenseMatrix64F gradient) {
        double[] actualization = new double[dim];
        cdi.getBranchActualization(getMatrixBufferOffsetIndex(nodeIndex), actualization);
        DenseMatrix64F qi = wrapDiagonal(actualization, 0, dim);

        DenseMatrix64F tmp = new DenseMatrix64F(dim, dim);
        DenseMatrix64F res = new DenseMatrix64F(dim, dim);
        CommonOps.mult(qi, gradient, tmp);
        CommonOps.multTransB(-1, tmp, qi, res);

        CommonOps.addEquals(res, gradient);

        double[] attenuation = elasticModel.getEigenValuesStrengthOfSelection();
        DenseMatrix64F variance = wrap(
                ((MultivariateIntegrator) cdi).getVariance(getEigenBufferOffsetIndex(0)),
                0, dim, dim);

        return chainRuleStationaryWrtAttenuationDiagonal(attenuation, variance, tmp, res);
    }

    private DenseMatrix64F chainRuleStationaryWrtAttenuationDiagonal(double[] attenuation, DenseMatrix64F variance,
                                                                     DenseMatrix64F tmp, DenseMatrix64F grad) {

        DenseMatrix64F res = new DenseMatrix64F(dim, 1);

        CommonOps.elementMult(variance, grad, tmp);

        for (int i = 0; i < dim; i++) {
            double sum = 0.0;
            for (int j = 0; j < dim; j++) {
                sum -= tmp.unsafe_get(i, j) * computeAttenuationFactor(attenuation, i, j);
            }
            res.unsafe_set(i, 0, sum);
        }

        return res;
    }

    private double computeAttenuationFactor(double[] attenuation, int i, int j) {
        return 2.0 / Math.pow(attenuation[i] + attenuation[j], 2);
    }

    DenseMatrix64F getGradientDisplacementWrtAttenuation(NodeRef node,
                                                         ContinuousDiffusionIntegrator cdi,
                                                         BranchSufficientStatistics statistics,
                                                         DenseMatrix64F gradient) {
        assert !tree.isRoot(node) : "Gradient wrt actualization is not available for the root.";
        if (hasDiagonalActualization()) {
            return getGradientDisplacementWrtAttenuationDiagonal(cdi, statistics, node.getNumber(), gradient);
        } else {
            throw new RuntimeException("not yet implemented");
        }
    }

    private DenseMatrix64F getGradientDisplacementWrtAttenuationDiagonal(ContinuousDiffusionIntegrator cdi, BranchSufficientStatistics statistics,
                                                                         int nodeIndex, DenseMatrix64F gradient) {
        // q_i
        double[] actualization = new double[dim];
        cdi.getBranchActualization(getMatrixBufferOffsetIndex(nodeIndex), actualization);
        DenseMatrix64F qi = wrapDiagonal(actualization, 0, dim);
        // q_i^-1
        double[] actualizationInv = new double[dim];
        for (int i = 0; i < dim; i++) {
            actualizationInv[i] = 1.0 / actualization[i];
        }
        DenseMatrix64F qiInv = wrapDiagonal(actualizationInv, 0, dim);
        // (I - q_i)^-1
        double[] actualizationIminusQInv = new double[dim];
        for (int i = 0; i < dim; i++) {
            actualizationIminusQInv[i] = 1.0 / (1.0 - actualization[i]);
        }
        DenseMatrix64F IminusqiInv = wrapDiagonal(actualizationIminusQInv, 0, dim);
        // ni
        DenseMatrix64F ni = statistics.getAbove().getRawMean();
        // beta_i
        double[] displacement = new double[dim];
        cdi.getBranchDisplacement(getMatrixBufferOffsetIndex(nodeIndex), displacement);
        DenseMatrix64F ri = wrap(displacement, 0, dim, 1);
        DenseMatrix64F betai = new DenseMatrix64F(dim, 1);
        CommonOps.mult(IminusqiInv, ri, betai);

        // factor
        DenseMatrix64F tmp = new DenseMatrix64F(dim, 1);
        DenseMatrix64F resFull = new DenseMatrix64F(dim, dim);
        DenseMatrix64F resDiag = new DenseMatrix64F(dim, 1);
        CommonOps.add(ni, -1, betai, resDiag);
        CommonOps.mult(qiInv, resDiag, tmp);
        CommonOps.multTransB(gradient, tmp, resFull);

        // Extract diag
        for (int i = 0; i < dim; ++i) {
            resDiag.unsafe_set(i, 0, resFull.unsafe_get(i, i));
        }

        // Wrt attenuation
        double ti = cdi.getBranchLength(getMatrixBufferOffsetIndex(nodeIndex));
        chainRuleActualizationWrtAttenuationDiagonal(qi, ti, resDiag);

        return resDiag;

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
                            var = Math.exp(-ep * ti) * Math.exp(-eq * tj) * ((Math.exp((ep + eq) * tij) - 1) / (ep + eq) + 1 / priorSampleSize) * traitVariance[p][q];
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
}