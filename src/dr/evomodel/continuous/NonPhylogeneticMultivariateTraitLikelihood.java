/*
 * NonPhylogeneticMultivariateTraitLikelihood.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.continuous;

import dr.evolution.tree.MutableTreeModel;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.math.distributions.WishartSufficientStatistics;
import dr.math.matrixAlgebra.Matrix;
import dr.math.matrixAlgebra.Vector;
import dr.util.Citable;

import java.util.List;
import java.util.logging.Logger;

/**
 * Integrated multivariate trait likelihood that assumes a fully-conjugate prior on the root and
 * no underlying tree structure.
 *
 * @author Gabriela Cybis
 * @author Marc A. Suchard
 * @author Bridgett vonHoldt
 */
public class NonPhylogeneticMultivariateTraitLikelihood extends FullyConjugateMultivariateTraitLikelihood {

    public NonPhylogeneticMultivariateTraitLikelihood(String traitName,
                                                     MutableTreeModel treeModel,
                                                     MultivariateDiffusionModel diffusionModel,
                                                     CompoundParameter traitParameter,
                                                     Parameter deltaParameter,
                                                     List<Integer> missingIndices,
                                                     boolean cacheBranches,
                                                     boolean scaleByTime,
                                                     boolean useTreeLength,
                                                     BranchRateModel rateModel,
                                                     Model samplingDensity,
                                                     boolean reportAsMultivariate,
                                                     double[] rootPriorMean,
                                                     double rootPriorSampleSize,
                                                      List<RestrictedPartials> partials,
                                                     boolean reciprocalRates,
                                                     boolean exchangeableTips) {
        super(traitName, treeModel, diffusionModel, traitParameter, deltaParameter, missingIndices, cacheBranches,
                scaleByTime, useTreeLength, rateModel, null, null, null, samplingDensity, reportAsMultivariate, rootPriorMean,
                partials, rootPriorSampleSize, reciprocalRates);
        this.exchangeableTips = exchangeableTips;
        this.zeroHeightTip = findZeroHeightTip(treeModel);
        printInformtion2();
    }

    private int findZeroHeightTip(Tree tree) {
        for (int i = 0; i < tree.getExternalNodeCount(); ++i) {
            NodeRef tip = tree.getExternalNode(i);
            if (tree.getNodeHeight(tip) == 0.0) {
                return i;
            }
        }
        return -1;
    }

    protected void printInformtion() {
        // Do nothing yet
    }

    protected void printInformtion2() {
        StringBuilder sb = new StringBuilder("Creating non-phylogenetic multivariate diffusion model:\n");
        sb.append("\tTrait: ").append(traitName).append("\n");
        sb.append("\tDiffusion process: ").append(diffusionModel.getId()).append("\n");
        sb.append("\tExchangeable tips: ").append((exchangeableTips ? "yes" : "no"));
        if (exchangeableTips) {
            sb.append(" initial inverse-weight = ").append(1.0 / getLengthToRoot(treeModel.getExternalNode(0)));
        }
        sb.append("\n");
        sb.append(extraInfo());
        sb.append("\tPlease cite:\n");
        sb.append(Citable.Utils.getCitationString(this));


        sb.append("\n\tDiffusion dimension   : ").append(dimTrait).append("\n");
        sb.append(  "\tNumber of observations: ").append(numData).append("\n");
        Logger.getLogger("dr.evomodel").info(sb.toString());
    }

    protected double getTreeLength() {
        double treeLength = 0;

        double rootHeight = treeModel.getNodeHeight(treeModel.getRoot());
        treeLength = 0;
        for (int i = 0; i < treeModel.getExternalNodeCount(); i++) {
            NodeRef node = treeModel.getExternalNode(i);
            treeLength += rootHeight - treeModel.getNodeHeight(node); // Bug was here
        }
        return treeLength;
    }
        
    private class SufficientStatistics {
        double sumWeight;
        double productWeight;
        double innerProduct;
        int nonMissingTips;

        SufficientStatistics(double sumWeight, double productWeight, double innerProduct, int nonMissingTips) {
            this.sumWeight = sumWeight;
            this.productWeight = productWeight;
            this.innerProduct = innerProduct;
            this.nonMissingTips = nonMissingTips;
        }
    }

    protected double getLengthToRoot(NodeRef nodeRef) {
        final double height;
        if (exchangeableTips) {
            height = getRescaledLengthToRoot(treeModel.getExternalNode(zeroHeightTip));
        } else {
            height = getRescaledLengthToRoot(nodeRef);
        }
        return height;
    }

    // Useful identity for computing outerproducts for Wishart statistics
    // \sum (y_i - \bar{y}) (y_i - \bar{y})^{t} = \sum y_i y_i^{t} - n \bar{y} \bar{y}^t

    private SufficientStatistics computeInnerProductsForTips(double[][] traitPrecision, double[] tmpVector) {

        // Compute the contribution of each datum at the root
        final int rootIndex = treeModel.getRoot().getNumber();
        final int meanOffset = dim * rootIndex;

        // Zero-out root mean
        for (int d = 0; d < dim; ++d) {
            meanCache[meanOffset + d] = 0;
        }

        double innerProducts = 0.0;

        // Compute the contribution of each datum at the root
        double productWeight = 1.0;
        double sumWeight = 0.0;
        int nonMissingTips = 0;
        for (int i = 0; i < treeModel.getExternalNodeCount(); ++i) {
            NodeRef tipNode = treeModel.getExternalNode(i);
            final int tipNumber = tipNode.getNumber();

            double tipWeight = 0.0;
            if (!missingTraits.isCompletelyMissing(tipNumber)) {
                
                tipWeight = 1.0 / getLengthToRoot(tipNode);

                int tipOffset = dim * tipNumber;
                int rootOffset = dim * rootIndex;

                for (int datum = 0; datum < numData; ++datum) {
                    // TODO Make faster when dimTrait == 1

                    // Add weighted tip value
                    for (int d = 0; d < dimTrait; ++d) {
                        meanCache[rootOffset + d] += tipWeight * meanCache[tipOffset + d];
                        tmpVector[d] = meanCache[tipOffset + d];
                    }

                    // Compute outer product
                    double yAy = computeWeightedAverageAndSumOfSquares(tmpVector, Ay, traitPrecision, dimTrait,
                                    tipWeight);
                    innerProducts += yAy;

                    if (DEBUG_NO_TREE) {
                        System.err.println("OP for " + tipNumber + " = " + yAy);
                        System.err.println("Value  = " + new Vector(tmpVector));
                        System.err.print  ("Prec   =\n" + new Matrix(traitPrecision));
                        System.err.println("weight = " + tipWeight + "\n");
                    }

                    tipOffset += dimTrait;
                    rootOffset += dimTrait;
                }

                if (computeWishartStatistics) {
                    incrementOuterProducts(tipNumber, tipWeight);
                }
            }

            if (tipWeight > 0.0) {
                sumWeight += tipWeight;
                productWeight *= tipWeight;
                ++nonMissingTips;
            }
        }

        lowerPrecisionCache[rootIndex] = sumWeight;
        normalize(meanCache, meanOffset, dim, sumWeight);

        if (computeWishartStatistics) {
            incrementOuterProducts(rootIndex, -sumWeight);
            wishartStatistics.incrementDf(-1);
        }

        return new SufficientStatistics(sumWeight, productWeight, innerProducts,
                nonMissingTips);
    }

    private void normalize(double[] x, int offset, int dim, double weight) {
        for (int d = 0; d < dim; ++d) {
            x[offset + d] /= weight;
        }
    }

    private void incrementOuterProducts(int nodeNumber, double nodeWeight) {

        final double[] outerProduct = wishartStatistics.getScaleMatrix();

        int tipOffset = dim * nodeNumber;
        for (int datum = 0; datum < numData; ++datum) {
            for (int i = 0; i < dim; ++i) {
                double yi = meanCache[tipOffset + i];
                for (int j = 0; j < dim; ++j) {
                    outerProduct[i * dim + j] += yi * meanCache[tipOffset +j] * nodeWeight;
                }
            }
            tipOffset += dimTrait;
        }

        wishartStatistics.incrementDf(1); // Peeled one node
    }

    protected boolean peel() {
        return false;
    }

    public double calculateLogLikelihood() {
        
        double[][] traitPrecision = diffusionModel.getPrecisionmatrix();
        double logDetTraitPrecision = Math.log(diffusionModel.getDeterminantPrecisionMatrix());
        double[] marginalRoot = tmp2;

        if (computeWishartStatistics) {
                wishartStatistics = new WishartSufficientStatistics(dimTrait);
        }

        // Compute the contribution of each datum at the root
        SufficientStatistics stats = computeInnerProductsForTips(traitPrecision, tmp2);

        double conditionalSumWeight = stats.sumWeight;
        double conditionalProductWeight = stats.productWeight;
        double innerProducts = stats.innerProduct;
        int nonMissingTips = stats.nonMissingTips;

        // Add in prior and integrate
        double sumWeight = conditionalSumWeight + rootPriorSampleSize;
        double productWeight = conditionalProductWeight * rootPriorSampleSize;
        double rootPrecision = productWeight / sumWeight;

        final int rootIndex = treeModel.getRoot().getNumber();
        int rootOffset = dim * rootIndex;

        for (int datum = 0; datum < numData; ++datum) {

            // Determine marginal root (scaled) mean
            for (int d = 0; d < dimTrait; ++d) {
                marginalRoot[d] = conditionalSumWeight * meanCache[rootOffset + d] + rootPriorSampleSize * rootPriorMean[d];
            }

            // Compute outer product contribution from prior
            double yAy1 = computeWeightedAverageAndSumOfSquares(rootPriorMean, Ay, traitPrecision, dimTrait,
                    rootPriorSampleSize);
            innerProducts += yAy1;  // TODO Only need to compute once

            if (DEBUG_NO_TREE) {
                System.err.println("OP for root");
                System.err.println("Value  = " + new Vector(rootPriorMean));
                System.err.print  ("Prec   = \n" + new Matrix(traitPrecision));
                System.err.println("Weight = " + rootPriorSampleSize + "\n");
            }

            // Compute outer product differences to complete square
            double yAy2 = computeWeightedAverageAndSumOfSquares(marginalRoot, Ay, traitPrecision, dimTrait,
                    1.0 / sumWeight);
            innerProducts -= yAy2;

            // Add prior on root contribution
            if (computeWishartStatistics) {

                final double[] outerProducts = wishartStatistics.getScaleMatrix();
                final double weight = conditionalSumWeight * rootPriorSampleSize / sumWeight;

                for (int i = 0; i < dimTrait; i++) {
                    final double diffi = meanCache[rootOffset + i] - rootPriorMean[i];
                    for (int j = 0; j < dimTrait; j++) {
                        outerProducts[i * dimTrait + j] += diffi * weight * (meanCache[rootOffset + j] - rootPriorMean[j]);
                    }
                }
                wishartStatistics.incrementDf(1);
            }

            rootOffset += dimTrait;
        }

        if (DEBUG_NO_TREE) {
            System.err.println("SumWeight    : " + sumWeight);
            System.err.println("ProductWeight: " + productWeight);
            System.err.println("Total OP     : " + innerProducts);
        }

        // Compute log likelihood
        double logLikelihood =
                -LOG_SQRT_2_PI * dimTrait * nonMissingTips * numData
                + 0.5 * logDetTraitPrecision * nonMissingTips * numData
                + 0.5 * Math.log(rootPrecision) * dimTrait * numData        
                - 0.5 * innerProducts;

        if (DEBUG_NO_TREE) {
            System.err.println("logLikelihood (final) = " + logLikelihood);
            System.err.println("numData = " + numData);
        }

        areStatesRedrawn = false;  // Should redraw internal node states when needed
        return logLikelihood;
    }

    private final boolean exchangeableTips;
    private final int zeroHeightTip;

    private static final boolean DEBUG_NO_TREE = false;
    private static final boolean NO_RESCALING = false;
}