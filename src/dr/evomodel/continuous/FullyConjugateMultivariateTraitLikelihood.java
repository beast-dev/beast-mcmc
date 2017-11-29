/*
 * FullyConjugateMultivariateTraitLikelihood.java
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
import dr.inference.model.*;
import dr.math.KroneckerOperation;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.math.distributions.WishartSufficientStatistics;
import dr.math.interfaces.ConjugateWishartStatisticsProvider;
import dr.math.matrixAlgebra.IllegalDimension;
import dr.math.matrixAlgebra.Matrix;
import dr.math.matrixAlgebra.Vector;
import dr.xml.Reportable;

import java.util.*;

/**
 * Integrated multivariate trait likelihood that assumes a fully-conjugate prior on the root.
 * The fully-conjugate prior is a multivariate normal distribution with a precision scaled by
 * diffusion process
 *
 * @author Marc A. Suchard
 */
public class FullyConjugateMultivariateTraitLikelihood extends IntegratedMultivariateTraitLikelihood
        implements ConjugateWishartStatisticsProvider, GibbsSampleFromTreeInterface, Reportable {

//    public FullyConjugateMultivariateTraitLikelihood(String traitName,
//                                                     MutableTreeModel treeModel,
//                                                     MultivariateDiffusionModel diffusionModel,
//                                                     CompoundParameter traitParameter,
//                                                     Parameter deltaParameter,
//                                                     List<Integer> missingIndices,
//                                                     boolean cacheBranches,
//                                                     boolean scaleByTime,
//                                                     boolean useTreeLength,
//                                                     BranchRateModel rateModel,
//                                                     Model samplingDensity,
//                                                     boolean reportAsMultivariate,
//                                                     double[] rootPriorMean,
//                                                     double rootPriorSampleSize,
//                                                     boolean reciprocalRates) {
//
//        super(traitName, treeModel, diffusionModel, traitParameter, deltaParameter, missingIndices, cacheBranches, scaleByTime,
//                useTreeLength, rateModel, samplingDensity, reportAsMultivariate, reciprocalRates);
//
//        // fully-conjugate multivariate normal with own mean and prior sample size
//        this.rootPriorMean = rootPriorMean;
//        this.rootPriorSampleSize = rootPriorSampleSize;
//
//        priorInformationKnown = false;
//    }
//
//
//    public FullyConjugateMultivariateTraitLikelihood(String traitName,
//                                                     MutableTreeModel treeModel,
//                                                     MultivariateDiffusionModel diffusionModel,
//                                                     CompoundParameter traitParameter,
//                                                     Parameter deltaParameter,
//                                                     List<Integer> missingIndices,
//                                                     boolean cacheBranches,
//                                                     boolean scaleByTime,
//                                                     boolean useTreeLength,
//                                                     BranchRateModel rateModel,
//                                                     List<BranchRateModel> driftModels,
//                                                     Model samplingDensity,
//                                                     boolean reportAsMultivariate,
//                                                     double[] rootPriorMean,
//                                                     double rootPriorSampleSize,
//                                                     boolean reciprocalRates) {
//
//        super(traitName, treeModel, diffusionModel, traitParameter, deltaParameter, missingIndices, cacheBranches, scaleByTime,
//                useTreeLength, rateModel, driftModels, samplingDensity, reportAsMultivariate, reciprocalRates);
//
//        // fully-conjugate multivariate normal with own mean and prior sample size
//        this.rootPriorMean = rootPriorMean;
//        this.rootPriorSampleSize = rootPriorSampleSize;
//
//        priorInformationKnown = false;
//    }

    public FullyConjugateMultivariateTraitLikelihood(String traitName,
                                                     MutableTreeModel treeModel,
                                                     MultivariateDiffusionModel diffusionModel,
                                                     CompoundParameter traitParameter,
                                                     Parameter deltaParameter,
                                                     List<Integer> missingIndices,
                                                     boolean cacheBranches,
                                                     boolean scaleByTime,
                                                     boolean useTreeLength,
                                                     BranchRateModel rateModel,
                                                     List<BranchRateModel> driftModels,
                                                     List<BranchRateModel> optimalValues,
                                                     BranchRateModel strengthOfSelection,
                                                     Model samplingDensity,
                                                     boolean reportAsMultivariate,
                                                     double[] rootPriorMean,
                                                     List<RestrictedPartials> partials,
                                                     double rootPriorSampleSize,
                                                     boolean reciprocalRates) {

        super(traitName, treeModel, diffusionModel, traitParameter, deltaParameter, missingIndices, cacheBranches, scaleByTime,
                useTreeLength, rateModel, driftModels, optimalValues, strengthOfSelection, samplingDensity, partials, reportAsMultivariate, reciprocalRates);

        // fully-conjugate multivariate normal with own mean and prior sample size
        this.rootPriorMean = rootPriorMean;
        this.rootPriorSampleSize = rootPriorSampleSize;

        priorInformationKnown = false;
    }


    //TODO temporary function so everything will compile. Need to actually write this.
    public FullyConjugateMultivariateTraitLikelihood semiClone(CompoundParameter traitParameter){
        return this;
    }

//    public FullyConjugateMultivariateTraitLikelihood semiClone(CompoundParameter traitParameter){
//        return new FullyConjugateMultivariateTraitLikelihood(this.traitName, this.treeModel, this.diffusionModel, traitParameter,
//                this.deltaParameter, this.missingIndices, this.cacheBranches, this.scaleByTime, this.useTreeLength, this.getBranchRateModel(),
//                this.optimalValues, this.strengthOfSelection, this.samplingDensity, this.reportAsMultivariate, this.rootPriorMean,
//                this.rootPriorSampleSize, this.reciprocalRates);
//    }

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

    @Override
    public MatrixParameterInterface getPrecisionParameter() {
        return diffusionModel.getPrecisionParameter();
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

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type){
        if(variable==traitParameter
                &&(Parameter.ChangeType.ADDED==type || Parameter.ChangeType.REMOVED==type)
                ){
            dimKnown = false;
            dim = traitParameter.getParameter(0).getDimension();
            numData = dim / getDimTrait();
            meanCache = new double[dim * treeModel.getNodeCount()];
//            storedMeanCache = new double[meanCache.length];
            drawnStates = new double[dim * treeModel.getNodeCount()];
        }
        PostPreKnown=false;
        super.handleVariableChangedEvent(variable,index,type);
    }

    @Override
    public void storeState() {
        storedNumData = numData;
        storedDim = dim;
        super.storeState();
        storedPostPreKnown=PostPreKnown;
        storedDimKnown=dimKnown;
        if(preP!=null)
         System.arraycopy(preP, 0, storedPreP, 0,preP.length);
        if(preMeans!=null){
            for(int i = 0; i < preMeans.length; i++)
                storedPreMeans[i] = preMeans[i].clone();
        }

    }

    @Override
    public void restoreState() {
        dim = storedDim;
        numData = storedNumData;
        drawnStates = new double[dim * treeModel.getNodeCount()];
        super.restoreState();
        PostPreKnown=storedPostPreKnown;
        priorInformationKnown = false;
        double[] tempPreP = storedPreP;
        storedPreP = preP;
        preP = tempPreP;
        preMeans=storedPreMeans;
        double[][] preMeansTemp = preMeans;
        preMeans = storedPreMeans;
        storedPreMeans = preMeansTemp;
        dimKnown=storedDimKnown;
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

    public void doPreOrderTraversal(NodeRef node) {

        if(preP==null){
            preP=new double[treeModel.getNodeCount()];
            storedPreP=new double[treeModel.getNodeCount()];
        }
        if(!dimKnown){
            preMeans=new double[treeModel.getNodeCount()][dim];
            storedPreMeans=new double[treeModel.getNodeCount()][dim];
            dimKnown=true;
        }

        final int thisNumber = node.getNumber();


        if (treeModel.isRoot(node)) {
            preP[thisNumber] = rootPriorSampleSize;
            for (int j = 0; j < dim; j++) {
                preMeans[thisNumber][j]
                        = rootPriorMean[j % dimTrait];
            }


        } else {

            final NodeRef parentNode = treeModel.getParent(node);
            final NodeRef sibNode = getSisterNode(node);

            final int parentNumber = parentNode.getNumber();
            final int sibNumber = sibNode.getNumber();



	/*

			  if (treeModel.isRoot(parentNode)){
				  //partial precisions
				    final double precisionParent = rootPriorSampleSize;
			        final double precisionSib = postP[sibNumber];
			        final double thisPrecision=1/treeModel.getBranchLength(node);
			        double tp= precisionParent + precisionSib;
			        preP[thisNumber]= tp*thisPrecision/(tp+thisPrecision);

			        //partial means

			        for (int j =0; j<dim;j++){
			        	preMeans[thisNumber][j] = (precisionParent*preMeans[parentNumber][j] + precisionSib*rootPriorMean[j])/(precisionParent+precisionSib);
			        }

			  }else{
	*/
            //partial precisions
            final double precisionParent = preP[parentNumber];
            final double precisionSib = upperPrecisionCache[sibNumber];
            final double thisPrecision = 1 / getRescaledBranchLengthForPrecision(node);
            double tp = precisionParent + precisionSib;
            preP[thisNumber] = tp * thisPrecision / (tp + thisPrecision);

            //partial means

            for (int j = 0; j < dim; j++) {
                preMeans[thisNumber][j] = (precisionParent * preMeans[parentNumber][j] + precisionSib * cacheHelper.getMeanCache()[sibNumber*dim+j]) / (precisionParent + precisionSib);
            }
        }

        if (treeModel.isExternal(node)) {
            return;
        } else {
            doPreOrderTraversal(treeModel.getChild(node, 0));
            doPreOrderTraversal(treeModel.getChild(node, 1));

        }

    }

    public NodeRef getSisterNode(NodeRef node) {
        NodeRef sib0 = treeModel.getChild(treeModel.getParent(node), 0);
        NodeRef sib1 = treeModel.getChild(treeModel.getParent(node), 1);


        if (sib0 == node) {
            return sib1;
        } else return sib0;

    }

    public double[] getConditionalMean(int taxa){
        setup();


//            double[] answer=new double[getRootNodeTrait().length];

        double[] mean = new double[dim];
        for (int i = 0; i < dim; i++) {
            mean[i] = preMeans[taxa][i];
        }

        return mean;
    }

    public double[][] getConditionalMeans(){
        setup();


//            double[] answer=new double[getRootNodeTrait().length];


        return preMeans;
    }

    public double getPrecisionFactor(int taxa){
        setup();
        return preP[taxa];
    }

    public double[] getPrecisionFactors(){
        setup();
        return preP;
    }

    public double[][] getConditionalPrecision(int taxa){
         setup();




        double[][] precisionParam =diffusionModel.getPrecisionmatrix();
//        double[][] answer=new double[getRootNodeTrait().length][ getRootNodeTrait().length];
        double p = getPrecisionFactor(taxa);

        double[][] thisP = new double[dim][dim];
        for (int i = 0; i < getNumData(); i++) {
            for (int j = 0; j < getDimTrait(); j++) {
                for (int k = 0; k < getDimTrait(); k++) {
//                System.out.println("P: "+p);
//                System.out.println("I: "+i+", J: "+j+" value:"+precisionParam[i][j]);
                    thisP[i * getDimTrait() + j][i * getDimTrait() + k] = p * precisionParam[j][k];

                }
            }
        }

        return thisP;

    }

    private void setup(){
        if(!PostPreKnown){
            double[][] traitPrecision = diffusionModel.getPrecisionmatrix();
            double logDetTraitPrecision = Math.log(diffusionModel.getDeterminantPrecisionMatrix());

            final boolean computeWishartStatistics = getComputeWishartSufficientStatistics();

            if (computeWishartStatistics) {
                wishartStatistics = new WishartSufficientStatistics(dimTrait);
            }

            // Use dynamic programming to compute conditional likelihoods at each internal node
            postOrderTraverse(treeModel, treeModel.getRoot(), traitPrecision, logDetTraitPrecision, computeWishartStatistics);

            doPreOrderTraversal(treeModel.getRoot());}
        PostPreKnown=true;

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

                final double[] outerProducts = wishartStatistics.getScaleMatrix();

                final double weight = conditionalRootPrecision * rootPriorSampleSize * marginalVariance;
                for (int i = 0; i < dimTrait; i++) {
                    final double diffi = conditionalRootMean[i] - rootPriorMean[i];
                    for (int j = 0; j < dimTrait; j++) {
                        outerProducts[i * dimTrait + j] += diffi * weight * (conditionalRootMean[j] - rootPriorMean[j]);
                    }
                }
                wishartStatistics.incrementDf(1);
            }
        } else {
            // 1D is very simple
            final double x = conditionalRootMean[0] * conditionalRootPrecision + rootPriorMean[0] * rootPriorSampleSize;
            square = x * x * treePrecisionMatrix[0][0] * marginalVariance;

            if (computeWishartStatistics) {
                final double[] outerProducts = wishartStatistics.getScaleMatrix();
                final double y = conditionalRootMean[0] - rootPriorMean[0];
                outerProducts[0] += y * y * conditionalRootPrecision * rootPriorSampleSize * marginalVariance;
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

    double[] preP;
    double[][] preMeans;

    double[] storedPreP;
    double[][] storedPreMeans;

    Boolean PostPreKnown=false;
    Boolean storedPostPreKnown=false;

    private boolean priorInformationKnown = false;
    private double zBz; // Prior sum-of-squares contribution

    private boolean dimKnown=false;
    private boolean storedDimKnown=false;
    private int storedDim;
    private int storedNumData;

    protected boolean computeWishartStatistics = false;
    private double[] ascertainedData = null;
    private static final boolean DEBUG_ASCERTAINMENT = false;

    private double vectorMin(double[] vec) {
        double min = Double.MAX_VALUE;
        for (int i = 0; i < vec.length; ++i) {
            min = Math.min(min, vec[i]);
        }
        return min;
    }

    private double matrixMin(double[][] mat) {
        double min = Double.MAX_VALUE;
        for (int i = 0; i < mat.length; ++i) {
            min = Math.min(min, vectorMin(mat[i]));
        }
        return min;
    }

    private double vectorMax(double[] vec) {
        double max = - Double.MAX_VALUE;
        for (int i = 0; i < vec.length; ++i) {
            max = Math.max(max, vec[i]);
        }
        return max;
    }

    private double matrixMax(double[][] mat) {
        double max = -Double.MAX_VALUE;
        for (int i = 0; i < mat.length; ++i) {
            max = Math.max(max, vectorMax(mat[i]));
        }
        return max;
    }

    private double vectorSum(double[] vec) {
        double sum = 0.0;
        for (int i = 0; i < vec.length; ++i) {
            sum += vec[i];
        }
        return sum;
    }

    private double matrixSum(double[][] mat) {
        double sum = 0.0;
        for (int i = 0; i < mat.length; ++i) {
            sum += vectorSum(mat[i]);
        }
        return sum;
    }

    @Override
    public String getReport() {
        StringBuilder sb = new StringBuilder();
//        sb.append(this.g)
//        System.err.println("Hello");
        sb.append("Tree:\n");
        sb.append(getId()).append("\t");
        sb.append(treeModel.toString());
        sb.append("\n\n");

        double[][] treeVariance = computeTreeVariance(true);
        double[][] traitPrecision = getDiffusionModel().getPrecisionmatrix();
        Matrix traitVariance = new Matrix(traitPrecision).inverse();

        double[][] jointVariance = KroneckerOperation.product(treeVariance, traitVariance.toComponents());

        sb.append("Tree variance:\n");
        sb.append(new Matrix(treeVariance));
        sb.append(matrixMin(treeVariance)).append("\t").append(matrixMax(treeVariance)).append("\t").append(matrixSum(treeVariance));
        sb.append("\n\n");
        sb.append("Trait variance:\n");
        sb.append(traitVariance);
        sb.append("\n\n");
//        sb.append("Joint variance:\n");
//        sb.append(new Matrix(jointVariance));
//        sb.append("\n\n");

        sb.append("Tree dim: " + treeVariance.length + "\n");
        sb.append("data dim: " + jointVariance.length);
        sb.append("\n\n");

        double[] data = new double[jointVariance.length];
        System.arraycopy(meanCache, 0, data, 0, jointVariance.length);

        if (nodeToClampMap != null) {
            int offset = treeModel.getExternalNodeCount() * getDimTrait();
            for(Map.Entry<NodeRef, RestrictedPartials> clamps : nodeToClampMap.entrySet()) {
                double[] partials = clamps.getValue().getPartials();
                for (int i = 0; i < partials.length; ++i) {
                    data[offset] = partials[i];
                    ++offset;
                }
            }
        }

        sb.append("Data:\n");
        sb.append(new Vector(data)).append("\n");
        sb.append(data.length).append("\t").append(vectorMin(data)).append("\t").append(vectorMax(data)).append("\t").append(vectorSum(data));
        sb.append(treeModel.getNodeTaxon(treeModel.getExternalNode(0)).getId());
        sb.append("\n\n");

        MultivariateNormalDistribution mvn = new MultivariateNormalDistribution(new double[data.length], new Matrix(jointVariance).inverse().toComponents());
        double logDensity = mvn.logPdf(data);
        sb.append("logLikelihood: " + getLogLikelihood() + " == " + logDensity + "\n\n");

        final WishartSufficientStatistics sufficientStatistics = getWishartStatistics();
        final double[] outerProducts = sufficientStatistics.getScaleMatrix();

        sb.append("Outer-products (DP):\n");
        sb.append(new Vector(outerProducts));
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


    class NodeToRootDistance {
        NodeRef node;
        double distance;

        NodeToRootDistance(NodeRef node, double distance) {
            this.node = node;
            this.distance = distance;
        }
    }

    class NodeToRootDistanceList extends ArrayList<NodeToRootDistance> {

        NodeToRootDistanceList(NodeToRootDistanceList parentList) {
            super(parentList);
        }

        NodeToRootDistanceList() {
            super();
        }
    }

    private void addNodeToList(final NodeRef thisNode, NodeToRootDistanceList parentList, NodeToRootDistanceList[] tipLists) {

        if (!treeModel.isRoot(thisNode)) {
            double increment = getRescaledBranchLengthForPrecision(thisNode);
            if (parentList.size() > 0) {
                increment += parentList.get(parentList.size() - 1).distance;
            }
            parentList.add(new NodeToRootDistance(thisNode, increment));
        }

        if (treeModel.isExternal(thisNode)) {
            tipLists[thisNode.getNumber()] =  parentList;
        } else { // recurse
            NodeToRootDistanceList shallowCopy = new NodeToRootDistanceList(parentList);
            addNodeToList(treeModel.getChild(thisNode, 0), shallowCopy, tipLists);
            addNodeToList(treeModel.getChild(thisNode, 1), parentList, tipLists);
        }
    }

    private double getTimeBetweenNodeToRootLists(List<NodeToRootDistance> x, List<NodeToRootDistance> y) {
        if (x.get(0) != y.get(0)) {
            return 0.0;
        }

        int index = 1;
        while (x.get(index) == y.get(index)) {
            ++index;
        }
        return x.get(index - 1).distance;
    }

    public double[][] computeTreeVariance2(boolean includeRoot) {

        final int tipCount = treeModel.getExternalNodeCount();
        double[][] variance = new double[tipCount][tipCount];

        NodeToRootDistanceList[] tipToRootDistances = new NodeToRootDistanceList[tipCount];

        // Recurse down tree to generate lists
        addNodeToList(treeModel.getRoot(), new NodeToRootDistanceList(), tipToRootDistances);

        for (int i = 0; i < tipCount; ++i) {
            // Fill in diagonal
            List<NodeToRootDistance> iList = tipToRootDistances[i];
            double marginalTime = iList.get(iList.size() - 1).distance;
            variance[i][i] = marginalTime;

            for (int j = i + 1; j < tipCount; ++j) {
                List<NodeToRootDistance> jList = tipToRootDistances[j];

                double time = getTimeBetweenNodeToRootLists(iList, jList);
                variance[j][i] = variance[i][j] = time;
            }
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

    public double[][] computeTreeVariance(boolean includeRoot) {
        final int tipCount = treeModel.getExternalNodeCount();
        int length = tipCount;

        boolean DO_CLAMP = true;
        if (DO_CLAMP && nodeToClampMap != null) {
            length += nodeToClampMap.size();
        }
//        System.exit(-1);

        double[][] variance = new double[length][length];

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

        if (DO_CLAMP && nodeToClampMap != null) {
            List<RestrictedPartials> partialsList = new ArrayList<RestrictedPartials>();
            for (Map.Entry<NodeRef, RestrictedPartials> keySet : nodeToClampMap.entrySet()) {
                partialsList.add(keySet.getValue());
            }

            for (int i = 0; i < partialsList.size(); ++i) {
                RestrictedPartials partials = partialsList.get(i);
                NodeRef node = partials.getNode();

                variance[tipCount + i][tipCount + i] = getRescaledLengthToRoot(node) +
                        1.0 / partials.getPriorSampleSize();

                for (int j = 0; j < tipCount; ++j) {
                    NodeRef friend = treeModel.getExternalNode(j);
                    NodeRef mrca = TreeUtils.getCommonAncestor(treeModel, node, friend);
                    variance[j][tipCount + i] = getRescaledLengthToRoot(mrca);

                }

                for (int j = 0; j < i; ++j) {
                    NodeRef friend = partialsList.get(j).getNode();
                    NodeRef mrca = TreeUtils.getCommonAncestor(treeModel, node, friend);
                    variance[tipCount + j][tipCount + i] = getRescaledLengthToRoot(mrca);
                }
            }
        }

        // Make symmetric
        for (int i = 0; i < length; i++) {
            for (int j = i + 1; j < length; j++) {
                variance[j][i] = variance[i][j];
            }
        }

//        if (DEBUG) {
//            System.err.println("");
//            System.err.println("New tree conditional variance:\n" + new Matrix(variance));
//        }
//
//        variance = removeMissingTipsInTreeVariance(variance); // Automatically prune missing tips
//
//        if (DEBUG) {
//            System.err.println("");
//            System.err.println("New tree (trimmed) conditional variance:\n" + new Matrix(variance));
//        }

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
        return TreeUtils.getCommonAncestorNode(treeModel, leafNames);
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
