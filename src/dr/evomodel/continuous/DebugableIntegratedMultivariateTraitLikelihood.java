/*
 * DebugableIntegratedMultivariateTraitLikelihood.java
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
import dr.evolution.tree.TreeUtils;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Model;
import dr.math.KroneckerOperation;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.math.matrixAlgebra.Matrix;
import dr.math.matrixAlgebra.SymmetricMatrix;
import dr.math.matrixAlgebra.Vector;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * This class contains all of the O(number of tips^2) functions to debug the dynamic programming in its superclasses
 * 
 * @author Marc A Suchard
 */
public class DebugableIntegratedMultivariateTraitLikelihood extends SemiConjugateMultivariateTraitLikelihood {

    public DebugableIntegratedMultivariateTraitLikelihood(String traitName,
                                                          MutableTreeModel treeModel,
                                                          MultivariateDiffusionModel diffusionModel,
                                                          CompoundParameter traitParameter,
                                                          List<Integer> missingIndices,
                                                          boolean cacheBranches,
                                                          boolean scaleByTime,
                                                          boolean useTreeLength,
                                                          BranchRateModel rateModel,
                                                          Model samplingDensity,
                                                          boolean reportAsMultivariate,
                                                          MultivariateNormalDistribution rootPrior,
                                                          boolean reciprocalRates,
                                                          List<RestrictedPartials> partials) {

        super(traitName, treeModel, diffusionModel, traitParameter, missingIndices, cacheBranches,
                scaleByTime, useTreeLength, rateModel, samplingDensity, reportAsMultivariate,
                rootPrior, reciprocalRates, partials);
    }

    protected double[] fillLeafTraits(int datum) {

        final int tipCount = treeModel.getExternalNodeCount();

        final int nonMissingTipCount = countNonMissingTips();

        double[] traits = new double[dimTrait * nonMissingTipCount];
        int index = 0;
        for (int i = 0; i < tipCount; i++) {

            if (!missingTraits.isCompletelyMissing(i)) {
                for (int k = 0; k < dimTrait; k++) {
                    traits[index++] = meanCache[dim * i + datum * dimTrait + k];
                }
            }
        }
        return traits;
    }

    protected double[][] removeMissingTipsInTreeVariance(double[][] variance) {

        final int tipCount = treeModel.getExternalNodeCount();
        final int nonMissing = countNonMissingTips();

        if (nonMissing == tipCount) { // Do nothing
            return variance;
        }

        double[][] outVariance = new double[nonMissing][nonMissing];

        int iReal = 0;
        for (int i = 0; i < tipCount; i++) {
            if (!missingTraits.isCompletelyMissing(i)) {

                int jReal = 0;
                for (int j = 0; j < tipCount; j++) {
                    if (!missingTraits.isCompletelyMissing(i)) {

                        outVariance[iReal][jReal] = variance[i][j];

                        jReal++;
                    }
                }
                iReal++;
            }
        }
        return outVariance;
    }

    protected double[][] computeTreeTraitPrecision(double[][] traitPrecision) {
        double[][] treePrecision = computeTreePrecision();
        if (dimTrait > 1) {
            treePrecision = KroneckerOperation.product(treePrecision, traitPrecision);
        } else {
            final double precision = traitPrecision[0][0];
            for (int i = 0; i < treePrecision.length; i++) {
                for (int j = 0; j < treePrecision[i].length; j++) {
                    treePrecision[i][j] *= precision;
                }
            }
        }
        return treePrecision;
    }

    public double[][] computeTreePrecision() {
        return new SymmetricMatrix(computeTreeVariance()).inverse().toComponents();
    }

    private NodeRef findMRCA(int iTip, int jTip) {
        Set<String> leafNames = new HashSet<String>();
        leafNames.add(treeModel.getTaxonId(iTip));
        leafNames.add(treeModel.getTaxonId(jTip));
        return TreeUtils.getCommonAncestorNode(treeModel, leafNames);
    }

    public int getNumberOfDatum() {
        return numData * countNonMissingTips();
    }

    protected double integrateLogLikelihoodAtRootFromFullTreeMatrix(double[][] treeTraitPrecisionMatrix,
                                                                    double[] tipTraits) {

        double logLikelihood = 0;
        final int tipCount = countNonMissingTips();

        // 1^t\Sigma^{-1} y + Pz
        double[] mean = Ay;
        for (int i = 0; i < dimTrait; i++) {
            mean[i] = 0;
            for (int j = 0; j < dimTrait; j++) {
                mean[i] += rootPriorPrecision[i][j] * rootPriorMean[j];
            }

            for (int j = 0; j < tipCount; j++) {
                final int rowOffset = j * dimTrait + i;
                for (int k = 0; k < tipCount * dimTrait; k++) {
                    mean[i] += treeTraitPrecisionMatrix[rowOffset][k] * tipTraits[k];
                }
            }
        }

        // 1^t \Sigma^{-1} 1 + P
        double[][] precision = tmpM;
        for (int i = 0; i < dimTrait; i++) {
            for (int j = 0; j < dimTrait; j++) {
                precision[i][j] = rootPriorPrecision[i][j];
                for (int k = 0; k < tipCount; k++) {
                    for (int l = 0; l < tipCount; l++) {
                        precision[i][j] += treeTraitPrecisionMatrix[k * dimTrait + i][l * dimTrait + j];
                    }
                }
            }
        }
        double[] normalizedMean = tmp2;
        double[][] variance = new SymmetricMatrix(precision).inverse().toComponents();
        for (int i = 0; i < dimTrait; i++) {
            normalizedMean[i] = 0.0;
            for (int j = 0; j < dimTrait; j++) {
                normalizedMean[i] += variance[i][j] * mean[j];
            }
        }
        mean = normalizedMean;

        // We know:  y ~ MVN(x, A) and x ~ N(m, B)
        // Therefore p(x | y) = N( (A+B)^{-1}(Ay + Bm), A + B)
        // We want: p( y ) = p( y | x ) p( x ) / p( x | y ) for any value x, say x = 0

        logLikelihood += MultivariateNormalDistribution.logPdf(
                rootPriorMean, new double[rootPriorMean.length], rootPriorPrecision,
                logRootPriorPrecisionDeterminant, 1.0
        );

        logLikelihood -= MultivariateNormalDistribution.logPdf(
                mean, new double[mean.length], precision,
                Math.log(MultivariateNormalDistribution.calculatePrecisionMatrixDeterminate(precision)), 1.0
        );

        if (DEBUG) {
            System.err.println("Mean = " + new Vector(mean));
            System.err.println("Prec = " + new Matrix(precision));
            System.err.println("log density = " + logLikelihood);
        }
        return logLikelihood;
    }

    public double[][] computeTreeVariance() {
        final int tipCount = treeModel.getExternalNodeCount();
        double[][] variance = new double[tipCount][tipCount];

        for (int i = 0; i < tipCount; i++) {

            // Fill in diagonal
            double marginalTime = getRescaledLengthToRoot(treeModel.getExternalNode(i));
            variance[i][i] = marginalTime;

            // Fill in upper right triangle,

            for (int j = i + 1; j < tipCount; j++) {
                NodeRef mrca = findMRCA(i, j);
                variance[i][j] = getRescaledLengthToRoot(mrca);
            }
        }

        // Make symmetric
        for (int i = 0; i < tipCount; i++) {
            for (int j = i + 1; j < tipCount; j++) {
                variance[j][i] = variance[i][j];
            }
        }

        if (DEBUG) {
            System.err.println("");
            System.err.println("New tree conditional variance:\n" + new Matrix(variance));
        }

        variance = removeMissingTipsInTreeVariance(variance); // Automatically prune missing tips

        if (DEBUG) {
            System.err.println("");
            System.err.println("New tree (trimmed) conditional variance:\n" + new Matrix(variance));
        }

        return variance;
    }

    protected int countNonMissingTips() {
        int tipCount = treeModel.getExternalNodeCount();
        for (int i = 0; i < tipCount; i++) {
            if (missingTraits.isCompletelyMissing(i)) {
                tipCount--;
            }
        }
        return tipCount;
    }

    public void checkViaLargeMatrixInversion() {

        // Perform a check based on filling in the (dimTrait * tipCount) * (dimTrait * tipCount) precision matrix
        // And then integrating out the root trait value

        // Form \Sigma^{-1} (precision) = (tree precision) %x% (trait precision)

        double[][] treeTraitPrecisionMatrix = computeTreeTraitPrecision(diffusionModel.getPrecisionmatrix());

        double totalLogDensity = 0;

        for (int datum = 0; datum < numData; datum++) {

            double[] tipTraits = fillLeafTraits(datum);

            System.err.println("Datum #" + datum);
            System.err.println("tipTraits = " + new Vector(tipTraits));
            System.err.println("tipPrecision = \n" + new Matrix(treeTraitPrecisionMatrix));

            double checkLogLikelihood = MultivariateNormalDistribution.logPdf(tipTraits, new double[tipTraits.length], treeTraitPrecisionMatrix,
                    Math.log(MultivariateNormalDistribution.calculatePrecisionMatrixDeterminate(treeTraitPrecisionMatrix)), 1.0);

            System.err.println("tipDensity = " + checkLogLikelihood + " (should match final likelihood when root not integrated out and no missing data)");

            // Convolve root prior
            if (integrateRoot) {
                checkLogLikelihood += integrateLogLikelihoodAtRootFromFullTreeMatrix(treeTraitPrecisionMatrix, tipTraits);
            }
            totalLogDensity += checkLogLikelihood;
        }
        System.err.println("Total logLikelihood (via tree) = " + totalLogDensity);
    }

    private double[][] computeTipTraitOuterProduct(int tip0, int tip1) {
        double[][] outerProduct = new double[dimTrait][dimTrait];

        final int offset0 = dim * tip0;
        final int offset1 = dim * tip1;

        for (int i = 0; i < dimTrait; i++) {
            for (int j = 0; j < dimTrait; j++) {
                for (int k = 0; k < numData; k++) {
                    outerProduct[i][j] += meanCache[offset0 + k * dimTrait + i] * meanCache[offset1 + k * dimTrait + j];
                }
            }
        }
        return outerProduct;
    }

    private void computeAllTipTraitOuterProducts() {
        final int nTips = treeModel.getExternalNodeCount();

        if (tipTraitOuterProducts == null) {
            tipTraitOuterProducts = new double[nTips][nTips][][];
        }

        for (int i = 0; i < nTips; i++) {
            if (!missingTraits.isCompletelyMissing(i)) {
                tipTraitOuterProducts[i][i] = computeTipTraitOuterProduct(i, i);
                for (int j = i + 1; j < nTips; j++) {
                    if (!missingTraits.isCompletelyMissing(i)) {
                        tipTraitOuterProducts[j][i] = tipTraitOuterProducts[i][j] = computeTipTraitOuterProduct(i, j);
                    } else {
                        tipTraitOuterProducts[j][i] = tipTraitOuterProducts[i][j] = null;
                    }
                }
            } else {
                for (int j = 0; j < nTips; j++) {
                    tipTraitOuterProducts[i][j] = null;
                }
            }
        }
    }

    // Returns the outer product of the tip traits for taxon 0 and taxon 1,
    // or null if either taxon 0 or taxon 1 is missing

    public double[][] getTipTraitOuterProduct(int tip0, int tip1) {
        if (updateOuterProducts) {
            computeAllTipTraitOuterProducts();
            updateOuterProducts = false;
        }
        return tipTraitOuterProducts[tip0][tip1];
    }

    protected boolean updateOuterProducts = true;
    protected double[][][][] tipTraitOuterProducts = null;
}
