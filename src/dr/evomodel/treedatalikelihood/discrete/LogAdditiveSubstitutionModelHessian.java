/*
 * LogAdditiveSubstitutionModelHessian.java
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

import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.substmodel.LogRateSubstitutionModel;
import dr.util.Transform;

import java.util.Arrays;

final class LogAdditiveSubstitutionModelHessian {

    private static final double EMPIRICAL_FISHER_RIDGE = 1.0E-12;

    interface BranchScoreProvider {
        double[] getBranchScores(int nodeIndex);
    }

    private final AbstractLogAdditiveSubstitutionModelGradient gradient;
    private final TreeTrait branchDifferentialTraitProvider;
    private final BranchScoreProvider branchScoreProvider = new BranchScoreProvider() {
        @Override
        public double[] getBranchScores(int nodeIndex) {
            return computeBranchScores(nodeIndex);
        }
    };

    private double[] rawBranchDifferentials;
    private double[] generator;
    private double[] frequencies;
    private double[] branchScores;
    private double[] accumulatedBranchDifferentials;
    private double[] empiricalFisher;
    private double[][] empiricalFisherMatrix;
    private boolean normalize;
    private Transform transform;
    private double rateScalar;

    LogAdditiveSubstitutionModelHessian(AbstractLogAdditiveSubstitutionModelGradient gradient,
                                        TreeTrait branchDifferentialTraitProvider) {
        this.gradient = gradient;
        this.branchDifferentialTraitProvider = branchDifferentialTraitProvider;
    }

    double[] getDiagonalHessianLogDensity() {
        requireBranchDifferentialTraitProvider("Diagonal");
        setUpBranchScoreState();
        empiricalFisher = getDiagonalHessianLogDensity(gradient.tree, gradient.getDimension(),
                branchScoreProvider, empiricalFisher);
        return empiricalFisher;
    }

    double[][] getHessianLogDensity() {
        requireBranchDifferentialTraitProvider("Full");
        setUpBranchScoreState();
        empiricalFisherMatrix = getHessianLogDensity(gradient.tree, gradient.getDimension(),
                branchScoreProvider, empiricalFisherMatrix);
        return empiricalFisherMatrix;
    }

    static double[] getDiagonalHessianLogDensity(Tree tree,
                                                 int dimension,
                                                 BranchScoreProvider branchScoreProvider,
                                                 double[] empiricalFisher) {
        if (empiricalFisher == null || empiricalFisher.length != dimension) {
            empiricalFisher = new double[dimension];
        } else {
            Arrays.fill(empiricalFisher, 0.0);
        }

        accumulateEmpiricalFisher(tree, dimension, branchScoreProvider, empiricalFisher, null);

        for (int i = 0; i < dimension; ++i) {
            empiricalFisher[i] = -(empiricalFisher[i] + EMPIRICAL_FISHER_RIDGE);
        }

        return empiricalFisher;
    }

    static double[][] getHessianLogDensity(Tree tree,
                                           int dimension,
                                           BranchScoreProvider branchScoreProvider,
                                           double[][] empiricalFisherMatrix) {
        if (empiricalFisherMatrix == null || empiricalFisherMatrix.length != dimension) {
            empiricalFisherMatrix = new double[dimension][dimension];
        } else {
            for (int i = 0; i < dimension; ++i) {
                Arrays.fill(empiricalFisherMatrix[i], 0.0);
            }
        }

        accumulateEmpiricalFisher(tree, dimension, branchScoreProvider, null, empiricalFisherMatrix);

        for (int i = 0; i < dimension; ++i) {
            empiricalFisherMatrix[i][i] -= EMPIRICAL_FISHER_RIDGE;
        }

        return empiricalFisherMatrix;
    }

    private double[] computeBranchScores(int nodeIndex) {
        final double[] branchDifferentials = accumulateBranchDifferentials(rawBranchDifferentials, nodeIndex);
        final double normalizationConstant = gradient.preProcessNormalization(branchDifferentials, generator, normalize);

        for (int i = 0; i < gradient.getDimension(); ++i) {
            branchScores[i] = gradient.processSingleGradientDimension(i, branchDifferentials, generator, frequencies,
                    normalize, normalizationConstant, rateScalar, transform, gradient.scaleRatesByFrequencies);
        }

        return branchScores;
    }

    private static void accumulateEmpiricalFisher(Tree tree,
                                                  int dimension,
                                                  BranchScoreProvider branchScoreProvider,
                                                  double[] diagonal,
                                                  double[][] full) {
        for (int nodeIndex = 0; nodeIndex < tree.getNodeCount(); ++nodeIndex) {
            if (tree.isRoot(tree.getNode(nodeIndex))) {
                continue;
            }

            final double[] scores = branchScoreProvider.getBranchScores(nodeIndex);
            if (diagonal != null) {
                for (int i = 0; i < dimension; ++i) {
                    diagonal[i] += scores[i] * scores[i];
                }
            } else {
                for (int i = 0; i < dimension; ++i) {
                    final double scoreI = scores[i];
                    for (int j = 0; j < dimension; ++j) {
                        full[i][j] -= scoreI * scores[j];
                    }
                }
            }
        }
    }

    private void setUpBranchScoreState() {
        rawBranchDifferentials = (double[]) branchDifferentialTraitProvider.getTrait(gradient.tree, null);

        final int generatorDimension = gradient.stateCount * gradient.stateCount;
        if (generator == null || generator.length != generatorDimension) {
            generator = new double[generatorDimension];
        }
        gradient.substitutionModel.getInfinitesimalMatrix(generator);

        frequencies = gradient.substitutionModel.getFrequencyModel().getFrequencies();
        normalize = gradient.substitutionModel.getNormalization();
        transform = (gradient.substitutionModel instanceof LogRateSubstitutionModel) ?
                ((LogRateSubstitutionModel) gradient.substitutionModel).getTransform() : null;
        rateScalar = normalize && transform != null ? 1 / gradient.substitutionModel.setupMatrix() : 0.0;

        final int dimension = gradient.getDimension();
        if (branchScores == null || branchScores.length != dimension) {
            branchScores = new double[dimension];
        }
    }

    private double[] accumulateBranchDifferentials(double[] branchDifferentials, int nodeIndex) {
        final int length = gradient.stateCount * gradient.stateCount;
        final int firstModel = gradient.crossProductAccumulationMap.get(0);
        final int modelStride = length;
        final int nodeStride = gradient.branchModel.getSubstitutionModels().size() * modelStride;
        final int nodeOffset = nodeIndex * nodeStride;

        if (accumulatedBranchDifferentials == null || accumulatedBranchDifferentials.length != length) {
            accumulatedBranchDifferentials = new double[length];
        } else {
            Arrays.fill(accumulatedBranchDifferentials, 0.0);
        }

        System.arraycopy(branchDifferentials, nodeOffset + firstModel * modelStride,
                accumulatedBranchDifferentials, 0, length);

        for (int i = 1; i < gradient.crossProductAccumulationMap.size(); ++i) {
            final int nextModel = gradient.crossProductAccumulationMap.get(i);
            final int offset = nodeOffset + nextModel * modelStride;
            for (int j = 0; j < length; ++j) {
                accumulatedBranchDifferentials[j] += branchDifferentials[offset + j];
            }
        }

        return accumulatedBranchDifferentials;
    }

    private void requireBranchDifferentialTraitProvider(String kind) {
        if (branchDifferentialTraitProvider == null) {
            throw new RuntimeException(kind + " empirical-Fisher preconditioning requires branch-wise " +
                    "substitution-model differentials, but this likelihood delegate does not provide them.");
        }
    }
}
