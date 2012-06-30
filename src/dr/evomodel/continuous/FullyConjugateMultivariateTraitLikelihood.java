/*
 * FullyConjugateMultivariateTraitLikelihood.java
 *
 * Copyright (c) 2002-2012 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.evolution.tree.NodeRef;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.math.distributions.WishartSufficientStatistics;
import dr.math.interfaces.ConjugateWishartStatisticsProvider;
import dr.math.matrixAlgebra.Matrix;
import dr.math.matrixAlgebra.Vector;

import java.util.List;

/**
 * Integrated multivariate trait likelihood that assumes a fully-conjugate prior on the root.
 * The fully-conjugate prior is a multivariate normal distribution with a precision scaled by
 * diffusion process
 *
 * @author Marc A. Suchard
 */
public class FullyConjugateMultivariateTraitLikelihood extends IntegratedMultivariateTraitLikelihood implements ConjugateWishartStatisticsProvider {

    public FullyConjugateMultivariateTraitLikelihood(String traitName,
                                                     TreeModel treeModel,
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


    protected double getRescaledLengthToRoot(NodeRef nodeRef) {

        double length = 0;
        NodeRef parent = treeModel.getParent(nodeRef);

        if (!treeModel.isRoot(parent)) {
            length += getRescaledLengthToRoot(parent);
        }

        length += getRescaledBranchLength(nodeRef);

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

        System.err.println("Got here subclass: " + loglikelihood);
        System.err.println("logValue         : " + (logRemainders + logPdf));
        System.err.println("");
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
}
