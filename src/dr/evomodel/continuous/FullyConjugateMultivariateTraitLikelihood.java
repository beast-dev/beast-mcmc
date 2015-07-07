/*
 * FullyConjugateMultivariateTraitLikelihood.java
 *
 * Copyright (c) 2002-2013 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.evolution.tree.MultivariateTraitTree;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.math.KroneckerOperation;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.math.distributions.NormalDistribution;
import dr.math.distributions.WishartSufficientStatistics;
import dr.math.interfaces.ConjugateWishartStatisticsProvider;
import dr.math.matrixAlgebra.IllegalDimension;
import dr.math.matrixAlgebra.Matrix;
import dr.math.matrixAlgebra.Vector;
import dr.xml.Reportable;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Integrated multivariate trait likelihood that assumes a fully-conjugate prior on the root.
 * The fully-conjugate prior is a multivariate normal distribution with a precision scaled by
 * diffusion process
 *
 * @author Marc A. Suchard
 */
public class FullyConjugateMultivariateTraitLikelihood extends IntegratedMultivariateTraitLikelihood implements ConjugateWishartStatisticsProvider, Reportable {

    public FullyConjugateMultivariateTraitLikelihood(String traitName,
                                                     MultivariateTraitTree treeModel,
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
                                                     boolean reciprocalRates) {

        super(traitName, treeModel, diffusionModel, traitParameter, deltaParameter, missingIndices, cacheBranches, scaleByTime,
                useTreeLength, rateModel, samplingDensity, reportAsMultivariate, reciprocalRates);

        // fully-conjugate multivariate normal with own mean and prior sample size
        this.rootPriorMean = rootPriorMean;
        this.rootPriorSampleSize = rootPriorSampleSize;

        priorInformationKnown = false;
    }


    public FullyConjugateMultivariateTraitLikelihood(String traitName,
                                                     MultivariateTraitTree treeModel,
                                                     MultivariateDiffusionModel diffusionModel,
                                                     CompoundParameter traitParameter,
                                                     Parameter deltaParameter,
                                                     List<Integer> missingIndices,
                                                     boolean cacheBranches,
                                                     boolean scaleByTime,
                                                     boolean useTreeLength,
                                                     BranchRateModel rateModel,
                                                     List<BranchRateModel> driftModels,
                                                     Model samplingDensity,
                                                     boolean reportAsMultivariate,
                                                     double[] rootPriorMean,
                                                     double rootPriorSampleSize,
                                                     boolean reciprocalRates) {

        super(traitName, treeModel, diffusionModel, traitParameter, deltaParameter, missingIndices, cacheBranches, scaleByTime,
                useTreeLength, rateModel, driftModels, samplingDensity, reportAsMultivariate, reciprocalRates);

        // fully-conjugate multivariate normal with own mean and prior sample size
        this.rootPriorMean = rootPriorMean;
        this.rootPriorSampleSize = rootPriorSampleSize;

        priorInformationKnown = false;
    }

    public FullyConjugateMultivariateTraitLikelihood(String traitName,
                                                     MultivariateTraitTree treeModel,
                                                     MultivariateDiffusionModel diffusionModel,
                                                     CompoundParameter traitParameter,
                                                     Parameter deltaParameter,
                                                     List<Integer> missingIndices,
                                                     boolean cacheBranches,
                                                     boolean scaleByTime,
                                                     boolean useTreeLength,
                                                     BranchRateModel rateModel,
                                                     List<BranchRateModel> optimalValues,
                                                     BranchRateModel strengthOfSelection,
                                                     Model samplingDensity,
                                                     boolean reportAsMultivariate,
                                                     double[] rootPriorMean,
                                                     double rootPriorSampleSize,
                                                     boolean reciprocalRates) {

        super(traitName, treeModel, diffusionModel, traitParameter, deltaParameter, missingIndices, cacheBranches, scaleByTime,
                useTreeLength, rateModel, optimalValues, strengthOfSelection, samplingDensity, reportAsMultivariate, reciprocalRates);

        // fully-conjugate multivariate normal with own mean and prior sample size
        this.rootPriorMean = rootPriorMean;
        this.rootPriorSampleSize = rootPriorSampleSize;

        priorInformationKnown = false;
    }


    public double getRescaledLengthToRoot(NodeRef nodeRef) {

        double length = 0;

        if (!treeModel.isRoot(nodeRef)) {
            NodeRef parent = treeModel.getParent(nodeRef);
            length += getRescaledBranchLengthForPrecision(nodeRef) + getRescaledLengthToRoot(parent);
        }

        return length;
    }


    protected double calculateAscertainmentCorrection(int taxonIndex) {

        NodeRef tip = treeModel.getNode(taxonIndex);
        int nodeIndex = treeModel.getNode(taxonIndex).getNumber();

        if (ascertainedData == null) { // Assumes that ascertained data are fixed
            ascertainedData = new double[dimTrait];
        }

//        diffusionModel.diffusionPrecisionMatrixParameter.setParameterValue(0,2); // For debugging non-1 values
        double[][] traitPrecision = diffusionModel.getPrecisionmatrix();
        double logDetTraitPrecision = Math.log(diffusionModel.getDeterminantPrecisionMatrix());

        double lengthToRoot = getRescaledLengthToRoot(tip);
        double marginalPrecisionScalar = 1.0 / lengthToRoot + rootPriorSampleSize;
        double logLikelihood = 0;

        for (int datum = 0; datum < numData; ++datum) {

            // Get observed trait value
            System.arraycopy(meanCache, nodeIndex * dim + datum * dimTrait, ascertainedData, 0, dimTrait);

            if (DEBUG_ASCERTAINMENT) {
                System.err.println("Datum #" + datum);
                System.err.println("Value: " + new Vector(ascertainedData));
                System.err.println("Cond : " + lengthToRoot);
                System.err.println("MargV: " + 1.0 / marginalPrecisionScalar);
                System.err.println("MargP: " + marginalPrecisionScalar);
                System.err.println("diffusion prec: " + new Matrix(traitPrecision));
            }

            double SSE;
            if (dimTrait > 1) {
                throw new RuntimeException("Still need to implement multivariate ascertainment correction");
            } else {
                double precision = traitPrecision[0][0] * marginalPrecisionScalar;

                SSE = ascertainedData[0] * precision * ascertainedData[0];
            }

            double thisLogLikelihood = -LOG_SQRT_2_PI * dimTrait
                    + 0.5 * (logDetTraitPrecision + dimTrait * Math.log(marginalPrecisionScalar) - SSE);

            if (DEBUG_ASCERTAINMENT) {
                System.err.println("LogLik: " + thisLogLikelihood);
                dr.math.distributions.NormalDistribution normal = new dr.math.distributions.NormalDistribution(0,
                        Math.sqrt(1.0 / (traitPrecision[0][0] * marginalPrecisionScalar)));
                System.err.println("TTTLik: " + normal.logPdf(ascertainedData[0]));
                if (datum >= 10) {
                    System.exit(-1);
                }
            }
            logLikelihood += thisLogLikelihood;
        }
        return logLikelihood;
    }

//    public double getRootPriorSampleSize() {
//        return rootPriorSampleSize;
//    }

//    public double[] getRootPriorMean() {
//        double[] out = new double[rootPriorMean.length];
//        System.arraycopy(rootPriorMean, 0, out, 0, out.length);
//        return out;
//    }

    public WishartSufficientStatistics getWishartStatistics() {
        computeWishartStatistics = true;
        calculateLogLikelihood();
        computeWishartStatistics = false;
        return wishartStatistics;
    }

//    private double getLogPrecisionDetermination() {
//        return Math.log(diffusionModel.getDeterminantPrecisionMatrix()) + dimTrait * Math.log(rootPriorSampleSize);
//    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {

        if (model == diffusionModel) {
            priorInformationKnown = false;
        }
        super.handleModelChangedEvent(model, object, index);
    }

    public void restoreState() {
        super.restoreState();
        priorInformationKnown = false;
    }

    public void makeDirty() {
        super.makeDirty();
        priorInformationKnown = false;
    }

    public double getPriorSampleSize() {
        return rootPriorSampleSize;
    }

    public double[] getPriorMean() {
        return rootPriorMean;
    }

    @Override
    public boolean getComputeWishartSufficientStatistics() {
        return computeWishartStatistics;
    }

    protected void checkLogLikelihood(double loglikelihood, double logRemainders,
                                      double[] conditionalRootMean, double conditionalRootPrecision,
                                      double[][] traitPrecision) {

//        System.err.println("root cmean    : " + new Vector(conditionalRootMean));
//        System.err.println("root cprec    : " + conditionalRootPrecision);
//        System.err.println("diffusion prec: " + new Matrix(traitPrecision));
//
//        System.err.println("prior mean    : " + new Vector(rootPriorMean));
//        System.err.println("prior prec    : " + rootPriorSampleSize);

        double upperPrecision = conditionalRootPrecision * rootPriorSampleSize / (conditionalRootPrecision + rootPriorSampleSize);
//        System.err.println("root cprec    : " + upperPrecision);

        double[][] newPrec = new double[traitPrecision.length][traitPrecision.length];

        for (int i = 0; i < traitPrecision.length; ++i) {
            for (int j = 0; j < traitPrecision.length; ++j) {
                newPrec[i][j] = traitPrecision[i][j] * upperPrecision;
            }
        }
        MultivariateNormalDistribution mvn = new MultivariateNormalDistribution(rootPriorMean, newPrec);
        double logPdf = mvn.logPdf(conditionalRootMean);

        if (Math.abs(loglikelihood - logRemainders - logPdf) > 1E-3) {

            System.err.println("Got here subclass: " + loglikelihood);
            System.err.println("logValue         : " + (logRemainders + logPdf));
            System.err.println("logRemainder = " + logRemainders);
            System.err.println("");
        }
//        System.err.println("logRemainders    : " + logRemainders);
//        System.err.println("logPDF           : " + logPdf);
//        System.exit(-1);
    }

    protected double integrateLogLikelihoodAtRoot(double[] conditionalRootMean,
                                                  double[] marginalRootMean,
                                                  double[][] notUsed,
                                                  double[][] treePrecisionMatrix, double conditionalRootPrecision) {
        final double square;
        final double marginalPrecision = conditionalRootPrecision + rootPriorSampleSize;
        final double marginalVariance = 1.0 / marginalPrecision;

        // square : (Ay + Bz)' (A+B)^{-1} (Ay + Bz)

        // A = conditionalRootPrecision * treePrecisionMatrix
        // B = rootPriorSampleSize * treePrecisionMatrix

        if (dimTrait > 1) {

            computeWeightedAverage(conditionalRootMean, 0, conditionalRootPrecision, rootPriorMean, 0, rootPriorSampleSize,
                    marginalRootMean, 0, dimTrait);

            square = computeQuadraticProduct(marginalRootMean, treePrecisionMatrix, marginalRootMean, dimTrait)
                    * marginalPrecision;

            if (computeWishartStatistics) {

                final double[][] outerProducts = wishartStatistics.getScaleMatrix();

                final double weight = conditionalRootPrecision * rootPriorSampleSize * marginalVariance;
                for (int i = 0; i < dimTrait; i++) {
                    final double diffi = conditionalRootMean[i] - rootPriorMean[i];
                    for (int j = 0; j < dimTrait; j++) {
                        outerProducts[i][j] += diffi * weight * (conditionalRootMean[j] - rootPriorMean[j]);
                    }
                }
                wishartStatistics.incrementDf(1);
            }
        } else {
            // 1D is very simple
            final double x = conditionalRootMean[0] * conditionalRootPrecision + rootPriorMean[0] * rootPriorSampleSize;
            square = x * x * treePrecisionMatrix[0][0] * marginalVariance;

            if (computeWishartStatistics) {
                final double[][] outerProducts = wishartStatistics.getScaleMatrix();
                final double y = conditionalRootMean[0] - rootPriorMean[0];
                outerProducts[0][0] += y * y * conditionalRootPrecision * rootPriorSampleSize * marginalVariance;
                wishartStatistics.incrementDf(1);
            }
        }

        if (!priorInformationKnown) {
            setRootPriorSumOfSquares(treePrecisionMatrix);
        }

        final double retValue = 0.5 * (dimTrait * Math.log(rootPriorSampleSize * marginalVariance) - zBz + square);

        if (DEBUG) {
            System.err.println("(Ay+Bz)(A+B)^{-1}(Ay+Bz) = " + square);
            System.err.println("density = " + retValue);
            System.err.println("zBz = " + zBz);
        }

        return retValue;
    }

    private void setRootPriorSumOfSquares(double[][] treePrecisionMatrix) {

        zBz = computeQuadraticProduct(rootPriorMean, treePrecisionMatrix, rootPriorMean, dimTrait) * rootPriorSampleSize;
        priorInformationKnown = true;
    }

    protected double[][] computeMarginalRootMeanAndVariance(double[] conditionalRootMean,
                                                            double[][] notUsed,
                                                            double[][] treeVarianceMatrix,
                                                            double conditionalRootPrecision) {

        final double[][] outVariance = tmpM; // Use a temporary buffer, will stay valid for only a short while

        computeWeightedAverage(conditionalRootMean, 0, conditionalRootPrecision, rootPriorMean, 0, rootPriorSampleSize,
                conditionalRootMean, 0, dimTrait);

        final double totalVariance = 1.0 / (conditionalRootPrecision + rootPriorSampleSize);
        for (int i = 0; i < dimTrait; i++) {
            for (int j = 0; j < dimTrait; j++) {
                outVariance[i][j] = treeVarianceMatrix[i][j] * totalVariance;
            }
        }

        return outVariance;
    }

    protected double[] rootPriorMean;
    protected double rootPriorSampleSize;

    private boolean priorInformationKnown = false;
    private double zBz; // Prior sum-of-squares contribution

    protected boolean computeWishartStatistics = false;
    private double[] ascertainedData = null;
    private static final boolean DEBUG_ASCERTAINMENT = false;

    @Override
    public String getReport() {
        StringBuilder sb = new StringBuilder();
//        sb.append(this.g)
//        System.err.println("Hello");
        sb.append("Tree:\n");
        sb.append(treeModel.toString());
        sb.append("\n\n");

        double[][] treeVariance = computeTreeVariance(true);
        double[][] traitPrecision = getDiffusionModel().getPrecisionmatrix();
        Matrix traitVariance = new Matrix(traitPrecision).inverse();

        double[][] jointVariance = KroneckerOperation.product(treeVariance, traitVariance.toComponents());

        sb.append("Tree variance:\n");
        sb.append(new Matrix(treeVariance));
        sb.append("\n\n");
        sb.append("Trait variance:\n");
        sb.append(traitVariance);
        sb.append("\n\n");
        sb.append("Joint variance:\n");
        sb.append(new Matrix(jointVariance));
        sb.append("\n\n");

        double[] data = new double[jointVariance.length];
        System.arraycopy(meanCache, 0, data, 0, jointVariance.length);

        sb.append("Data:\n");
        sb.append(new Vector(data));
        sb.append("\n\n");

        MultivariateNormalDistribution mvn = new MultivariateNormalDistribution(new double[data.length], new Matrix(jointVariance).inverse().toComponents());
        double logDensity = mvn.logPdf(data);
        sb.append("logLikelihood: " + getLogLikelihood() + " == " + logDensity + "\n\n");

        final WishartSufficientStatistics sufficientStatistics = getWishartStatistics();
        final double[][] outerProducts = sufficientStatistics.getScaleMatrix();

        sb.append("Outer-products (DP):\n");
        sb.append(new Matrix(outerProducts));
        sb.append(sufficientStatistics.getDf() + "\n");

        Matrix treePrecision = new Matrix(treeVariance).inverse();
        final int n = data.length / traitPrecision.length;
        final int p = traitPrecision.length;
        double[][] tmp = new double[n][p];

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < p; ++j) {
                tmp[i][j] = data[i * p + j];
            }
        }
        Matrix y = new Matrix(tmp);

        Matrix S = null;
        try {
            S = y.transpose().product(treePrecision).product(y); // Using Matrix-Normal form
        } catch (IllegalDimension illegalDimension) {
            illegalDimension.printStackTrace();
        }
        sb.append("Outer-products (from tree variance:\n");
        sb.append(S);
        sb.append("\n\n");

        return sb.toString();
    }


    private double[][] computeTreeVariance(boolean includeRoot) {
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

        if (includeRoot) {
            for (int i = 0; i < variance.length; ++i) {
                for (int j = 0; j < variance[i].length; ++j) {
                    variance[i][j] += 1.0 / getPriorSampleSize();
                }
            }
        }

        return variance;
    }

    private NodeRef findMRCA(int iTip, int jTip) {
        Set<String> leafNames = new HashSet<String>();
        leafNames.add(treeModel.getTaxonId(iTip));
        leafNames.add(treeModel.getTaxonId(jTip));
        return Tree.Utils.getCommonAncestorNode(treeModel, leafNames);
    }

    private double[][] removeMissingTipsInTreeVariance(double[][] variance) {

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

    private int countNonMissingTips() {
        int tipCount = treeModel.getExternalNodeCount();
        for (int i = 0; i < tipCount; i++) {
            if (missingTraits.isCompletelyMissing(i)) {
                tipCount--;
            }
        }
        return tipCount;
    }


}
