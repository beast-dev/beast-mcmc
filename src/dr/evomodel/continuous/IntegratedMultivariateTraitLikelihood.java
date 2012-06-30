/*
 * IntegratedMultivariateTraitLikelihood.java
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
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.loggers.LogColumn;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.math.distributions.WishartSufficientStatistics;
import dr.math.matrixAlgebra.Matrix;
import dr.math.matrixAlgebra.SymmetricMatrix;
import dr.math.matrixAlgebra.Vector;
import dr.util.Author;
import dr.util.Citation;

import java.util.Arrays;
import java.util.List;
import java.util.logging.Logger;

/**
 * A multivariate trait likelihood that analytically integrates out the unobserved trait values at all internal
 * and root nodes
 *
 * @author Marc A. Suchard
 */
public abstract class IntegratedMultivariateTraitLikelihood extends AbstractMultivariateTraitLikelihood {

    public static final double LOG_SQRT_2_PI = 0.5 * Math.log(2 * Math.PI);

    public IntegratedMultivariateTraitLikelihood(String traitName,
                                                 TreeModel treeModel,
                                                 MultivariateDiffusionModel diffusionModel,
                                                 CompoundParameter traitParameter,
                                                 List<Integer> missingIndices,
                                                 boolean cacheBranches, boolean scaleByTime, boolean useTreeLength,
                                                 BranchRateModel rateModel, Model samplingDensity,
                                                 boolean reportAsMultivariate,
                                                 boolean reciprocalRates) {

        this(traitName, treeModel, diffusionModel, traitParameter, null, missingIndices, cacheBranches, scaleByTime,
                useTreeLength, rateModel, samplingDensity, reportAsMultivariate, reciprocalRates);
    }

    public IntegratedMultivariateTraitLikelihood(String traitName,
                                                 TreeModel treeModel,
                                                 MultivariateDiffusionModel diffusionModel,
                                                 CompoundParameter traitParameter,
                                                 Parameter deltaParameter,
                                                 List<Integer> missingIndices,
                                                 boolean cacheBranches, boolean scaleByTime, boolean useTreeLength,
                                                 BranchRateModel rateModel, Model samplingDensity,
                                                 boolean reportAsMultivariate,
                                                 boolean reciprocalRates) {

        super(traitName, treeModel, diffusionModel, traitParameter, deltaParameter, missingIndices, cacheBranches, scaleByTime,
                useTreeLength, rateModel, samplingDensity, reportAsMultivariate, reciprocalRates);


        meanCache = new double[dim * treeModel.getNodeCount()];
        drawnStates = new double[dim * treeModel.getNodeCount()];
        upperPrecisionCache = new double[treeModel.getNodeCount()];
        lowerPrecisionCache = new double[treeModel.getNodeCount()];
        logRemainderDensityCache = new double[treeModel.getNodeCount()];

        if (cacheBranches) {
            storedMeanCache = new double[dim * treeModel.getNodeCount()];
            storedUpperPrecisionCache = new double[treeModel.getNodeCount()];
            storedLowerPrecisionCache = new double[treeModel.getNodeCount()];
            storedLogRemainderDensityCache = new double[treeModel.getNodeCount()];
        }

        missing = new boolean[treeModel.getNodeCount()];
        Arrays.fill(missing, true); // All internal and root nodes are missing

        // Set up reusable temporary storage
        Ay = new double[dimTrait];
        tmpM = new double[dimTrait][dimTrait];
        tmp2 = new double[dimTrait];

        zeroDimVector = new double[dim];

        setTipDataValuesForAllNodes(missingIndices);

    }

    private void setTipDataValuesForAllNodes(List<Integer> missingIndices) {
        for (int i = 0; i < treeModel.getExternalNodeCount(); i++) {
            NodeRef node = treeModel.getExternalNode(i);
            setTipDataValuesForNode(node);
        }
        for (Integer i : missingIndices) {
            int whichTip = i / dim;
            Logger.getLogger("dr.evomodel").info(
                    "\tMarking taxon " + treeModel.getTaxonId(whichTip) + " as completely missing");
            missing[whichTip] = true;
        }
    }

    public double getTotalTreePrecision() {
        getLogLikelihood(); // Do peeling if necessary
        final int rootIndex = treeModel.getRoot().getNumber();
        return lowerPrecisionCache[rootIndex];
    }

    private void setTipDataValuesForNode(NodeRef node) {
        // Set tip data values
        int index = node.getNumber();
        double[] traitValue = traitParameter.getParameter(index).getParameterValues();
        if (traitValue.length < dim) {
            throw new RuntimeException("The trait parameter for the tip with index, " + index + ", is too short");
        }
        System.arraycopy(traitValue, 0, meanCache, dim * index, dim);
        missing[index] = false;
    }

    public double[] getTipDataValues(int index) {
        double[] traitValue = new double[dim];
        System.arraycopy(meanCache, dim * index, traitValue, 0, dim);
        return traitValue;
    }

    public void setTipDataValuesForNode(int index, double[] traitValue) {
        // Set tip data values
        System.arraycopy(traitValue, 0, meanCache, dim * index, dim);
        makeDirty();
    }

    protected String extraInfo() {
        return "\tSample internal node traits: false\n";
    }

    public List<Citation> getCitations() {
        List<Citation> citations = super.getCitations();
        citations.add(
                new Citation(
                        new Author[]{
                                new Author("O", "Pybus"),
                                new Author("P", "Lemey"),
                                new Author("A", "Rambaut"),
                                new Author("MA", "Suchard")
                        },
                        Citation.Status.IN_PREPARATION
                )
        );
        return citations;
    }

    public double getLogDataLikelihood() {
        return getLogLikelihood();
    }

    public abstract boolean getComputeWishartSufficientStatistics();

    public double calculateLogLikelihood() {

        double logLikelihood = 0;
        double[][] traitPrecision = diffusionModel.getPrecisionmatrix();
        double logDetTraitPrecision = Math.log(diffusionModel.getDeterminantPrecisionMatrix());
        double[] conditionalRootMean = tmp2;

        final boolean computeWishartStatistics = getComputeWishartSufficientStatistics();

        if (computeWishartStatistics) {
//            if (wishartStatistics == null) {
            wishartStatistics = new WishartSufficientStatistics(dimTrait);
//            } else {
//                wishartStatistics.clear();
//            }
        }

        // Use dynamic programming to compute conditional likelihoods at each internal node
        postOrderTraverse(treeModel, treeModel.getRoot(), traitPrecision, logDetTraitPrecision, computeWishartStatistics);

        if (DEBUG) {
            System.err.println("mean: " + new Vector(meanCache));
            System.err.println("upre: " + new Vector(upperPrecisionCache));
            System.err.println("lpre: " + new Vector(lowerPrecisionCache));
            System.err.println("cach: " + new Vector(logRemainderDensityCache));
        }

        // Compute the contribution of each datum at the root
        final int rootIndex = treeModel.getRoot().getNumber();

        // Precision scalar of datum conditional on root
        double conditionalRootPrecision = lowerPrecisionCache[rootIndex];

        for (int datum = 0; datum < numData; datum++) {

            double thisLogLikelihood = 0;

            // Get conditional mean of datum conditional on root
            System.arraycopy(meanCache, rootIndex * dim + datum * dimTrait, conditionalRootMean, 0, dimTrait);

            if (DEBUG) {
                System.err.println("Datum #" + datum);
                System.err.println("root mean: " + new Vector(conditionalRootMean));
                System.err.println("root prec: " + conditionalRootPrecision);
                System.err.println("diffusion prec: " + new Matrix(traitPrecision));
            }

            // B = root prior precision
            // z = root prior mean
            // A = likelihood precision
            // y = likelihood mean

            // y'Ay
            double yAy = computeWeightedAverageAndSumOfSquares(conditionalRootMean, Ay, traitPrecision, dimTrait,
                    conditionalRootPrecision); // Also fills in Ay

            if (conditionalRootPrecision != 0) {
                thisLogLikelihood += -LOG_SQRT_2_PI * dimTrait
                        + 0.5 * (logDetTraitPrecision + dimTrait * Math.log(conditionalRootPrecision) - yAy);
            }

            if (DEBUG) {
                double[][] T = new double[dimTrait][dimTrait];
                for (int i = 0; i < dimTrait; i++) {
                    for (int j = 0; j < dimTrait; j++) {
                        T[i][j] = traitPrecision[i][j] * conditionalRootPrecision;
                    }
                }
                System.err.println("Conditional root MVN precision = \n" + new Matrix(T));
                System.err.println("Conditional root MVN density = " + MultivariateNormalDistribution.logPdf(
                        conditionalRootMean, new double[dimTrait], T,
                        Math.log(MultivariateNormalDistribution.calculatePrecisionMatrixDeterminate(T)), 1.0));
            }

            if (integrateRoot) {
                // Integrate root trait out against rootPrior
                thisLogLikelihood += integrateLogLikelihoodAtRoot(conditionalRootMean, Ay, tmpM, traitPrecision,
                        conditionalRootPrecision); // Ay is destroyed
            }

            if (DEBUG) {
                System.err.println("yAy = " + yAy);
                System.err.println("logLikelihood (before remainders) = " + thisLogLikelihood +
                        " (should match conditional root MVN density when root not integrated out)");
            }

            logLikelihood += thisLogLikelihood;
        }

        logLikelihood += sumLogRemainders();

        if (DEBUG) { // Root trait is univariate!!!
            System.err.println("logLikelihood (final) = " + logLikelihood);
//            checkViaLargeMatrixInversion();
        }

        if (DEBUG_PNAS) {
            checkLogLikelihood(logLikelihood, sumLogRemainders(), conditionalRootMean,
                    conditionalRootPrecision, traitPrecision);
        }

        areStatesRedrawn = false;  // Should redraw internal node states when needed
        return logLikelihood;
    }

    protected void checkLogLikelihood(double loglikelihood, double logRemainders,
                                      double[] conditionalRootMean, double conditionalRootPrecision,
                                      double[][] traitPrecision) {
        // Do nothing; for checking PNAS paper
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if (variable == traitParameter) { // A tip value got updated
            if (index > dimTrait * treeModel.getExternalNodeCount()) {
                throw new RuntimeException("Attempting to update an invalid index");
            }
            meanCache[index] = traitParameter.getValue(index);
            likelihoodKnown = false;
        }
        super.handleVariableChangedEvent(variable, index, type);
    }

    protected static double computeWeightedAverageAndSumOfSquares(double[] y, double[] Ay, double[][] A,
                                                                  int dim, double scale) {
        // returns Ay and yAy
        double yAy = 0;
        for (int i = 0; i < dim; i++) {
            Ay[i] = 0;
            for (int j = 0; j < dim; j++)
                Ay[i] += A[i][j] * y[j] * scale;
            yAy += y[i] * Ay[i];
        }
        return yAy;
    }

    private double sumLogRemainders() {
        double sumLogRemainders = 0;
        for (double r : logRemainderDensityCache)
            sumLogRemainders += r;
        // Could skip leafs
        return sumLogRemainders;
    }

    protected abstract double integrateLogLikelihoodAtRoot(double[] conditionalRootMean,
                                                           double[] marginalRootMean,
                                                           double[][] temporaryStorage,
                                                           double[][] treePrecisionMatrix,
                                                           double conditionalRootPrecision);

    public void makeDirty() {
        super.makeDirty();
        areStatesRedrawn = false;
    }

    void postOrderTraverse(TreeModel treeModel, NodeRef node, double[][] precisionMatrix,
                           double logDetPrecisionMatrix, boolean cacheOuterProducts) {

        final int thisNumber = node.getNumber();

        if (treeModel.isExternal(node)) {

            // Fill in precision scalar, traitValues already filled in

            if (missing[thisNumber]) {
                upperPrecisionCache[thisNumber] = 0;
                lowerPrecisionCache[thisNumber] = 0; // Needed in the pre-order traversal
            } else { // not missing tip trait
                upperPrecisionCache[thisNumber] = 1.0 / getRescaledBranchLength(node);
                lowerPrecisionCache[thisNumber] = Double.POSITIVE_INFINITY;
            }
            return;
        }

        final NodeRef childNode0 = treeModel.getChild(node, 0);
        final NodeRef childNode1 = treeModel.getChild(node, 1);

        postOrderTraverse(treeModel, childNode0, precisionMatrix, logDetPrecisionMatrix, cacheOuterProducts);
        postOrderTraverse(treeModel, childNode1, precisionMatrix, logDetPrecisionMatrix, cacheOuterProducts);

        final int childNumber0 = childNode0.getNumber();
        final int childNumber1 = childNode1.getNumber();
        final int meanOffset0 = dim * childNumber0;
        final int meanOffset1 = dim * childNumber1;
        final int meanThisOffset = dim * thisNumber;

        final double precision0 = upperPrecisionCache[childNumber0];
        final double precision1 = upperPrecisionCache[childNumber1];
        final double totalPrecision = precision0 + precision1;

        lowerPrecisionCache[thisNumber] = totalPrecision;

        // Multiple child0 and child1 densities

        if (totalPrecision == 0) {
            System.arraycopy(zeroDimVector, 0, meanCache, meanThisOffset, dim);
        } else {

//            computeWeightedMeanCache(meanThisOffset, meanOffset0, meanOffset1, precision0, precision1);

            computeWeightedAverage(
                    meanCache, meanOffset0, precision0,
                    meanCache, meanOffset1, precision1,
                    meanCache, meanThisOffset, dim);
        }

        if (!treeModel.isRoot(node)) {
            // Integrate out trait value at this node
            double thisPrecision = 1.0 / getRescaledBranchLength(node);
            upperPrecisionCache[thisNumber] = totalPrecision * thisPrecision / (totalPrecision + thisPrecision);
        }

        // Compute logRemainderDensity

        logRemainderDensityCache[thisNumber] = 0;

        if (precision0 != 0 && precision1 != 0) {

            incrementRemainderDensities(
                    precisionMatrix,
                    logDetPrecisionMatrix, thisNumber, meanThisOffset,
                    meanOffset0,
                    meanOffset1,
                    precision0,
                    precision1,
                    cacheOuterProducts);
        }
    }

    private void incrementRemainderDensities(double[][] precisionMatrix,
                                             double logDetPrecisionMatrix,
                                             int thisIndex,
                                             int thisOffset,
                                             int childOffset0,
                                             int childOffset1,
                                             double precision0,
                                             double precision1,
                                             boolean cacheOuterProducts) {

        final double remainderPrecision = precision0 * precision1 / (precision0 + precision1);

        if (cacheOuterProducts) {
            incrementOuterProducts(thisOffset, childOffset0, childOffset1, precision0, precision1);
        }

        for (int k = 0; k < numData; k++) {

            double childSS0 = 0;
            double childSS1 = 0;
            double crossSS = 0;

            for (int i = 0; i < dimTrait; i++) {

                final double wChild0i = meanCache[childOffset0 + k * dimTrait + i] * precision0;
                final double wChild1i = meanCache[childOffset1 + k * dimTrait + i] * precision1;

                for (int j = 0; j < dimTrait; j++) {

                    final double child0j = meanCache[childOffset0 + k * dimTrait + j];
                    final double child1j = meanCache[childOffset1 + k * dimTrait + j];

                    childSS0 += wChild0i * precisionMatrix[i][j] * child0j;
                    childSS1 += wChild1i * precisionMatrix[i][j] * child1j;

                    crossSS += (wChild0i + wChild1i) * precisionMatrix[i][j] * meanCache[thisOffset + k * dimTrait + j];
                }
            }

            logRemainderDensityCache[thisIndex] +=
                    -dimTrait * LOG_SQRT_2_PI
                            + 0.5 * (dimTrait * Math.log(remainderPrecision) + logDetPrecisionMatrix)
                            - 0.5 * (childSS0 + childSS1 - crossSS);
        }
    }

    private void incrementOuterProducts(int thisOffset,
                                        int childOffset0,
                                        int childOffset1,
                                        double precision0,
                                        double precision1) {

        final double[][] outerProduct = wishartStatistics.getScaleMatrix();

        for (int k = 0; k < numData; k++) {

            for (int i = 0; i < dimTrait; i++) {

                final double wChild0i = meanCache[childOffset0 + k * dimTrait + i] * precision0;
                final double wChild1i = meanCache[childOffset1 + k * dimTrait + i] * precision1;

                for (int j = 0; j < dimTrait; j++) {

                    final double child0j = meanCache[childOffset0 + k * dimTrait + j];
                    final double child1j = meanCache[childOffset1 + k * dimTrait + j];

                    outerProduct[i][j] += wChild0i * child0j;
                    outerProduct[i][j] += wChild1i * child1j;

                    outerProduct[i][j] -= (wChild0i + wChild1i) * meanCache[thisOffset + k * dimTrait + j];
                }
            }
        }
        wishartStatistics.incrementDf(1); // Peeled one node
    }

//    private void computeWeightedMeanCache(int thisOffset,
//                                          int childOffset0,
//                                          int childOffset1,
//                                          double precision0,
//                                          double precision1) {
//
//        final double totalVariance = 1.0 / (precision0 + precision1);
//        for (int i = 0; i < dim; i++) {
//            meanCache[thisOffset + i] = (meanCache[childOffset0 + i] * precision0 +
//                    meanCache[childOffset1 + i] * precision1)
//                    * totalVariance;
//        }
//    }

    protected double[] getRootNodeTrait() {
        return getTraitForNode(treeModel, treeModel.getRoot(), traitName);
    }

    public double[] getTraitForNode(Tree tree, NodeRef node, String traitName) {

//        if (tree != treeModel) {
//            throw new RuntimeException("Can only reconstruct states on treeModel given to constructor");
//        }

        getLogLikelihood();

        if (!areStatesRedrawn)
            redrawAncestralStates();

        int index = node.getNumber();

        double[] trait = new double[dim];
        System.arraycopy(drawnStates, index * dim, trait, 0, dim);
        return trait;
    }

    public void redrawAncestralStates() {

        double[][] treePrecision = diffusionModel.getPrecisionmatrix();
        double[][] treeVariance = new SymmetricMatrix(treePrecision).inverse().toComponents();

        preOrderTraverseSample(treeModel, treeModel.getRoot(), 0, treePrecision, treeVariance);

        if (DEBUG) {
            System.err.println("all draws = " + new Vector(drawnStates));
        }

        areStatesRedrawn = true;
    }

    public void storeState() {
        super.storeState();

        if (cacheBranches) {
            System.arraycopy(meanCache, 0, storedMeanCache, 0, meanCache.length);
            System.arraycopy(upperPrecisionCache, 0, storedUpperPrecisionCache, 0, upperPrecisionCache.length);
            System.arraycopy(lowerPrecisionCache, 0, storedLowerPrecisionCache, 0, lowerPrecisionCache.length);
            System.arraycopy(logRemainderDensityCache, 0, storedLogRemainderDensityCache, 0, logRemainderDensityCache.length);
        }
    }

    public void restoreState() {
        super.restoreState();

        if (cacheBranches) {
            double[] tmp;

            tmp = storedMeanCache;
            storedMeanCache = meanCache;
            meanCache = tmp;

            tmp = storedUpperPrecisionCache;
            storedUpperPrecisionCache = upperPrecisionCache;
            upperPrecisionCache = tmp;

            tmp = storedLowerPrecisionCache;
            storedLowerPrecisionCache = lowerPrecisionCache;
            lowerPrecisionCache = tmp;

            tmp = storedLogRemainderDensityCache;
            storedLogRemainderDensityCache = logRemainderDensityCache;
            logRemainderDensityCache = tmp;
        }
    }


    // Computes x^t A y, used many times in these computations

    protected static double computeQuadraticProduct(double[] x, double[][] A, double[] y, int dim) {
        double sum = 0;
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                sum += x[i] * A[i][j] * y[j];
            }
        }
        return sum;
    }

    // Computes the weighted average of two vectors, used many times in these computations

    protected static void computeWeightedAverage(double[] in0, int offset0, double weight0,
                                                 double[] in1, int offset1, double weight1,
                                                 double[] out2, int offset2,
                                                 int length) {

        final double totalInverseWeight = 1.0 / (weight0 + weight1);
        for (int i = 0; i < length; i++) {
            out2[offset2 + i] = (in0[offset0 + i] * weight0 + in1[offset1 + i] * weight1) * totalInverseWeight;
        }
    }

    protected abstract double[][] computeMarginalRootMeanAndVariance(double[] conditionalRootMean,
                                                                     double[][] treePrecisionMatrix,
                                                                     double[][] treeVarianceMatrix,
                                                                     double conditionalRootPrecision);


    private void preOrderTraverseSample(TreeModel treeModel, NodeRef node, int parentIndex, double[][] treePrecision,
                                        double[][] treeVariance) {

        final int thisIndex = node.getNumber();

        if (treeModel.isRoot(node)) {
            // draw root

            double[] rootMean = new double[dimTrait];
            final int rootIndex = treeModel.getRoot().getNumber();
            double rootPrecision = lowerPrecisionCache[rootIndex];

            for (int datum = 0; datum < numData; datum++) {
                System.arraycopy(meanCache, thisIndex * dim + datum * dimTrait, rootMean, 0, dimTrait);

                double[][] variance = computeMarginalRootMeanAndVariance(rootMean, treePrecision, treeVariance,
                        rootPrecision);

                double[] draw = MultivariateNormalDistribution.nextMultivariateNormalVariance(rootMean, variance);

                if (DEBUG_PREORDER) {
                    Arrays.fill(draw, 1.0);
                }

                System.arraycopy(draw, 0, drawnStates, rootIndex * dim + datum * dimTrait, dimTrait);

                if (DEBUG) {
                    System.err.println("Root mean: " + new Vector(rootMean));
                    System.err.println("Root var : " + new Matrix(variance));
                    System.err.println("Root draw: " + new Vector(draw));
                }
            }
        } else { // draw conditional on parentState

            if (!missing[thisIndex]) {

                System.arraycopy(meanCache, thisIndex * dim, drawnStates, thisIndex * dim, dim);

            } else {

                // This code should work for sampling a missing tip trait as well, but needs testing

                // parent trait at drawnStates[parentOffset]
                double precisionToParent = 1.0 / getRescaledBranchLength(node);
                double precisionOfNode = lowerPrecisionCache[thisIndex];
                double totalPrecision = precisionOfNode + precisionToParent;

                double[] mean = Ay; // temporary storage
                double[][] var = tmpM; // temporary storage

                for (int datum = 0; datum < numData; datum++) {

                    int parentOffset = parentIndex * dim + datum * dimTrait;
                    int thisOffset = thisIndex * dim + datum * dimTrait;

                    if (DEBUG) {
                        double[] parentValue = new double[dimTrait];
                        System.arraycopy(drawnStates, parentOffset, parentValue, 0, dimTrait);
                        System.err.println("Parent draw: " + new Vector(parentValue));
                        if (parentValue[0] != drawnStates[parentOffset]) {
                            throw new RuntimeException("Error in setting indices");
                        }
                    }

                    for (int i = 0; i < dimTrait; i++) {
                        mean[i] = (drawnStates[parentOffset + i] * precisionToParent
                                + meanCache[thisOffset + i] * precisionOfNode) / totalPrecision;
                        for (int j = 0; j < dimTrait; j++) {
                            var[i][j] = treeVariance[i][j] / totalPrecision;
                        }
                    }
                    double[] draw = MultivariateNormalDistribution.nextMultivariateNormalVariance(mean, var);
                    System.arraycopy(draw, 0, drawnStates, thisOffset, dimTrait);

                    if (DEBUG) {
                        System.err.println("Int prec: " + totalPrecision);
                        System.err.println("Int mean: " + new Vector(mean));
                        System.err.println("Int var : " + new Matrix(var));
                        System.err.println("Int draw: " + new Vector(draw));
                        System.err.println("");
                    }
                }
            }
        }

        if (peel() && !treeModel.isExternal(node)) {
            preOrderTraverseSample(treeModel, treeModel.getChild(node, 0), thisIndex, treePrecision, treeVariance);
            preOrderTraverseSample(treeModel, treeModel.getChild(node, 1), thisIndex, treePrecision, treeVariance);
        }
    }

    protected boolean peel() {
        return true;
    }

    public LogColumn[] getColumns() {
        return new LogColumn[]{
                new LikelihoodColumn(getId())};
    }

    protected boolean areStatesRedrawn = false;

    protected double[] meanCache;
    protected double[] upperPrecisionCache;
    protected double[] lowerPrecisionCache;
    private double[] logRemainderDensityCache;

    protected boolean[] missing;

    private double[] storedMeanCache;
    private double[] storedUpperPrecisionCache;
    private double[] storedLowerPrecisionCache;
    private double[] storedLogRemainderDensityCache;

    private double[] drawnStates;

    protected final boolean integrateRoot = true; // Set to false if conditioning on root value (not fully implemented)
    protected static boolean DEBUG = false;
    protected static boolean DEBUG_PREORDER = false;
    protected static boolean DEBUG_PNAS = false;

    private double[] zeroDimVector;

    protected WishartSufficientStatistics wishartStatistics;

    // Reusable temporary storage
    protected double[] Ay;
    protected double[][] tmpM;
    protected double[] tmp2;

}
