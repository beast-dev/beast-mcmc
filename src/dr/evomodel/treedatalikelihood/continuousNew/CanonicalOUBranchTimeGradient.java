/*
 * CanonicalOUBranchTimeGradient.java
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

package dr.evomodel.treedatalikelihood.continuousNew;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.DifferentiableBranchRates;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.RateRescalingScheme;
import dr.evomodel.treedatalikelihood.continuous.ContinuousRateTransformation;

import java.util.Arrays;

final class CanonicalOUBranchTimeGradient {

    private CanonicalOUBranchTimeGradient() {
    }

    static double[] branchRateParameterGradient(final Tree tree,
                                                final DifferentiableBranchRates branchRateModel,
                                                final ContinuousRateTransformation rateTransformation,
                                                final double[] effectiveTimeGradient) {
        final double[] branchRateGradient = new double[tree.getNodeCount() - 1];
        final double normalization = normalization(rateTransformation);

        for (int nodeIndex = 0; nodeIndex < tree.getNodeCount(); nodeIndex++) {
            final NodeRef node = tree.getNode(nodeIndex);
            if (tree.isRoot(node)) {
                continue;
            }

            final int destinationIndex = branchRateModel.getParameterIndexFromNode(node);
            if (destinationIndex < 0 || destinationIndex >= branchRateGradient.length) {
                throw new IllegalArgumentException(
                        "Branch-rate model mapped node " + node.getNumber()
                                + " to invalid branch-gradient index " + destinationIndex + ".");
            }

            branchRateGradient[destinationIndex] += effectiveTimeGradient[nodeIndex]
                    * tree.getBranchLength(node)
                    * normalization
                    * branchRateModel.getBranchRateDifferential(tree, node);
        }

        return branchRateModel.updateGradientLogDensity(
                branchRateGradient, null, 0, branchRateGradient.length);
    }

    static void fillRawBranchLengthGradient(final Tree tree,
                                            final BranchRateModel branchRateModel,
                                            final ContinuousRateTransformation rateTransformation,
                                            final double[] effectiveTimeGradient,
                                            final double[] rawBranchLengthGradient) {
        Arrays.fill(rawBranchLengthGradient, 0.0);
        final double normalization = normalization(rateTransformation);
        if (scheme(rateTransformation) == RateRescalingScheme.TREE_HEIGHT) {
            throw new UnsupportedOperationException(
                    "Raw branch-length gradients with tree-height normalization are parameterization-dependent; "
                            + "use node-height gradients instead.");
        }
        final double normalizationDerivative;
        if (scheme(rateTransformation) == RateRescalingScheme.TREE_LENGTH) {
            normalizationDerivative = -normalization * normalization;
        } else {
            normalizationDerivative = 0.0;
        }
        final double dLogLdNormalization = dLogLdNormalization(tree, branchRateModel, effectiveTimeGradient);

        for (int nodeIndex = 0; nodeIndex < tree.getNodeCount(); nodeIndex++) {
            final NodeRef node = tree.getNode(nodeIndex);
            if (tree.isRoot(node)) {
                continue;
            }
            rawBranchLengthGradient[nodeIndex] =
                    effectiveTimeGradient[nodeIndex] * branchRate(tree, branchRateModel, node) * normalization
                            + dLogLdNormalization * normalizationDerivative;
        }
    }

    static void fillNodeHeightGradient(final TreeModel tree,
                                       final BranchRateModel branchRateModel,
                                       final ContinuousRateTransformation rateTransformation,
                                       final double[] effectiveTimeGradient,
                                       final double[] nodeHeightGradient) {
        Arrays.fill(nodeHeightGradient, 0.0);

        final double normalization = normalization(rateTransformation);
        final RateRescalingScheme scheme = scheme(rateTransformation);
        final double dLogLdNormalization = dLogLdNormalization(tree, branchRateModel, effectiveTimeGradient);
        final double normalizationSquared = normalization * normalization;

        for (int internalIndex = 0; internalIndex < tree.getInternalNodeCount(); internalIndex++) {
            final NodeRef node = tree.getNode(internalIndex + tree.getExternalNodeCount());
            double gradient = 0.0;

            for (int childIndex = 0; childIndex < tree.getChildCount(node); childIndex++) {
                final NodeRef child = tree.getChild(node, childIndex);
                gradient += effectiveTimeGradient[child.getNumber()]
                        * branchRate(tree, branchRateModel, child)
                        * normalization;
            }

            if (!tree.isRoot(node)) {
                gradient -= effectiveTimeGradient[node.getNumber()]
                        * branchRate(tree, branchRateModel, node)
                        * normalization;
            }

            if (scheme == RateRescalingScheme.TREE_LENGTH) {
                final double dTreeLengthDh =
                        tree.getChildCount(node) - (tree.isRoot(node) ? 0.0 : 1.0);
                gradient -= dLogLdNormalization * normalizationSquared * dTreeLengthDh;
            } else if (scheme == RateRescalingScheme.TREE_HEIGHT && tree.isRoot(node)) {
                gradient -= dLogLdNormalization * normalizationSquared;
            }

            nodeHeightGradient[internalIndex] = gradient;
        }
    }

    static double branchRate(final Tree tree,
                             final BranchRateModel branchRateModel,
                             final NodeRef node) {
        return branchRateModel == null ? 1.0 : branchRateModel.getBranchRate(tree, node);
    }

    private static double dLogLdNormalization(final Tree tree,
                                              final BranchRateModel branchRateModel,
                                              final double[] effectiveTimeGradient) {
        double result = 0.0;
        for (int nodeIndex = 0; nodeIndex < tree.getNodeCount(); nodeIndex++) {
            final NodeRef node = tree.getNode(nodeIndex);
            if (!tree.isRoot(node)) {
                result += effectiveTimeGradient[nodeIndex]
                        * tree.getBranchLength(node)
                        * branchRate(tree, branchRateModel, node);
            }
        }
        return result;
    }

    private static double normalization(final ContinuousRateTransformation rateTransformation) {
        return rateTransformation == null ? 1.0 : rateTransformation.getNormalization();
    }

    private static RateRescalingScheme scheme(final ContinuousRateTransformation rateTransformation) {
        return rateTransformation == null
                ? RateRescalingScheme.NONE
                : rateTransformation.getRateRescalingScheme();
    }
}
