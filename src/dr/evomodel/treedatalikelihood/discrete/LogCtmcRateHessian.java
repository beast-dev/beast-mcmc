/*
 * LogCtmcRateHessian.java
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

package dr.evomodel.treedatalikelihood.discrete;

import dr.evolution.tree.TreeTrait;
import dr.evomodel.substmodel.GlmSubstitutionModel;
import dr.evomodel.substmodel.LogRateSubstitutionModel;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.GradientDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.inference.hmc.HessianWrtParameterProvider;

import java.util.Arrays;

public class LogCtmcRateHessian extends LogCtmcRateGradient implements HessianWrtParameterProvider {
// extending LogCtmcRateGradient to satisfy the gradient requirements in HessianWrtParameterProvider
    private final LogAdditiveSubstitutionModelHessian hessian;
    private final TreeTrait branchLogRateScoreTraitProvider;
    private final LogAdditiveSubstitutionModelHessian.BranchScoreProvider branchScoreProvider =
            new LogAdditiveSubstitutionModelHessian.BranchScoreProvider() {
                @Override
                public double[] getBranchScores(int nodeIndex) {
                    return accumulateBranchScores(rawBranchScores, nodeIndex, branchScoreDimension);
                }
            };

    private double[] accumulatedBranchScores;
    private double[] empiricalFisher;
    private double[][] empiricalFisherMatrix;
    private double[] rawBranchScores;
    private int branchScoreDimension;

    public LogCtmcRateHessian(String traitName,
                              TreeDataLikelihood treeDataLikelihood,
                              BeagleDataLikelihoodDelegate likelihoodDelegate,
                              GlmSubstitutionModel substitutionModel) {
        super(traitName, treeDataLikelihood, likelihoodDelegate, substitutionModel);
        this.hessian = new LogAdditiveSubstitutionModelHessian(this,
                treeDataLikelihood.getTreeTrait(SpectralExactGradientDelegate.getBranchDifferentialTraitName(traitName)));
        this.branchLogRateScoreTraitProvider =
                treeDataLikelihood.getTreeTrait(SpectralExactGradientDelegate.getBranchLogRateScoreTraitName(traitName));
    }

    public LogCtmcRateHessian(String traitName,
                              TreeDataLikelihood treeDataLikelihood,
                              GradientDataLikelihoodDelegate likelihoodDelegate,
                              LogRateSubstitutionModel substitutionModel) {
        this(traitName, treeDataLikelihood, likelihoodDelegate, substitutionModel,
                ApproximationMode.FIRST_ORDER, false);
    }

    public LogCtmcRateHessian(String traitName,
                              TreeDataLikelihood treeDataLikelihood,
                              GradientDataLikelihoodDelegate likelihoodDelegate,
                              LogRateSubstitutionModel substitutionModel,
                              ApproximationMode mode,
                              boolean forceAllReal) {
        super(traitName, treeDataLikelihood, likelihoodDelegate, substitutionModel, mode, forceAllReal);
        this.hessian = new LogAdditiveSubstitutionModelHessian(this,
                treeDataLikelihood.getTreeTrait(SpectralExactGradientDelegate.getBranchDifferentialTraitName(traitName)));
        this.branchLogRateScoreTraitProvider =
                treeDataLikelihood.getTreeTrait(SpectralExactGradientDelegate.getBranchLogRateScoreTraitName(traitName));
    }

    @Override
    public double[] getDiagonalHessianLogDensity() {

        if (branchLogRateScoreTraitProvider == null ||
                getParameter().getDimension() != stateCount * (stateCount - 1)) {
            return hessian.getDiagonalHessianLogDensity();
        }

        final int dimension = getParameter().getDimension();
        rawBranchScores = (double[]) branchLogRateScoreTraitProvider.getTrait(tree, null);
        branchScoreDimension = dimension;
        empiricalFisher = LogAdditiveSubstitutionModelHessian.getDiagonalHessianLogDensity(
                tree, dimension, branchScoreProvider, empiricalFisher);

        return empiricalFisher;
    }

    @Override
    public double[][] getHessianLogDensity() {

        if (branchLogRateScoreTraitProvider == null ||
                getParameter().getDimension() != stateCount * (stateCount - 1)) {
            return hessian.getHessianLogDensity();
        }

        final int dimension = getParameter().getDimension();
        rawBranchScores = (double[]) branchLogRateScoreTraitProvider.getTrait(tree, null);
        branchScoreDimension = dimension;
        empiricalFisherMatrix = LogAdditiveSubstitutionModelHessian.getHessianLogDensity(
                tree, dimension, branchScoreProvider, empiricalFisherMatrix);

        return empiricalFisherMatrix;
    }

    private double[] accumulateBranchScores(double[] branchScores, int nodeIndex, int dimension) {
        final int firstModel = crossProductAccumulationMap.get(0);
        final int modelStride = dimension;
        final int nodeStride = branchModel.getSubstitutionModels().size() * modelStride;
        final int nodeOffset = nodeIndex * nodeStride;

        if (accumulatedBranchScores == null || accumulatedBranchScores.length != dimension) {
            accumulatedBranchScores = new double[dimension];
        } else {
            Arrays.fill(accumulatedBranchScores, 0.0);
        }

        System.arraycopy(branchScores, nodeOffset + firstModel * modelStride,
                accumulatedBranchScores, 0, dimension);

        for (int i = 1; i < crossProductAccumulationMap.size(); ++i) {
            final int nextModel = crossProductAccumulationMap.get(i);
            final int offset = nodeOffset + nextModel * modelStride;
            for (int j = 0; j < dimension; ++j) {
                accumulatedBranchScores[j] += branchScores[offset + j];
            }
        }

        return accumulatedBranchScores;
    }
}
