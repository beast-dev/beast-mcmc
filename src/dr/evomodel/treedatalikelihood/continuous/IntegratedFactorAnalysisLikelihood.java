/*
 * IntegratedFactorTraitDataModel.java
 *
 * Copyright (c) 2002-2018 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.evolution.tree.*;
import dr.evomodel.treedatalikelihood.continuous.cdi.PrecisionType;
import dr.evomodel.treedatalikelihood.preorder.ContinuousExtensionDelegate;
import dr.evomodel.treedatalikelihood.preorder.ModelExtensionProvider;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.inference.model.*;
import dr.math.KroneckerOperation;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.math.matrixAlgebra.IllegalDimension;
import dr.math.matrixAlgebra.Matrix;
import dr.math.matrixAlgebra.Vector;
import dr.math.matrixAlgebra.WrappedVector;
import dr.math.matrixAlgebra.missingData.InversionResult;
import dr.math.matrixAlgebra.missingData.MissingOps;
import dr.util.TaskPool;
import dr.xml.*;
import org.ejml.data.DenseMatrix64F;

import java.util.*;

import static dr.evomodelxml.treelikelihood.TreeTraitParserUtilities.STANDARDIZE;
import static dr.evomodelxml.treelikelihood.TreeTraitParserUtilities.TARGET_SD;
import static dr.math.matrixAlgebra.missingData.MissingOps.*;

/**
 * @author Marc A. Suchard
 * @author Gabriel Hassler
 */

public class IntegratedFactorAnalysisLikelihood extends AbstractModelLikelihood
        implements ContinuousTraitPartialsProvider, ModelExtensionProvider.NormalExtensionProvider, Reportable {

    private final int[] fullyObservedTraits;
    private final int[] partiallyMissingTraits;
    private boolean observedInnerProductKnown = false;
    private final DenseMatrix64F observedInnerProduct;
    // TODO: caching observedInnerProduct

    private static final PrecisionType precisionType = PrecisionType.FULL;
    private boolean[] missingTraitIndicators = null;

    private String tipTraitName;


    public IntegratedFactorAnalysisLikelihood(String name,
                                              CompoundParameter traitParameter,
                                              boolean[] missingIndicators,
                                              MatrixParameterInterface loadings,
                                              Parameter traitPrecision,
                                              double nuggetPrecision,
                                              TaskPool taskPool,
                                              CacheProvider cacheProvider) {

        super(name);

        this.traitParameter = traitParameter;
        this.loadingsTransposed = loadings;
        this.traitPrecision = traitPrecision;

        this.numTaxa = traitParameter.getParameterCount();
        this.dimTrait = traitParameter.getParameter(0).getDimension();
        this.numFactors = loadings.getColumnDimension();

        assert (dimTrait == loadings.getRowDimension());

        this.dimPartial = precisionType.getPartialsDimension(numFactors);

        addVariable(traitParameter);
        addVariable(loadings);
        addVariable(traitPrecision);

        this.missingDataIndicator = missingIndicators;
        this.missingDataIndices = ContinuousTraitPartialsProvider.indicatorToIndices(missingIndicators); //TODO: deprecate
        this.observedIndicators = setupObservedIndicators(missingDataIndices, numTaxa, dimTrait);
        this.observedDimensions = setupObservedDimensions(observedIndicators);

        List<Integer> observedList = new ArrayList<>();
        List<Integer> partialList = new ArrayList<>();
        setupObservedTraits(observedList, partialList);
        this.fullyObservedTraits = new int[observedList.size()];
        for (int i = 0; i < observedList.size(); i++) fullyObservedTraits[i] = observedList.get(i);
        this.partiallyMissingTraits = new int[partialList.size()];
        for (int i = 0; i < partialList.size(); i++) partiallyMissingTraits[i] = partialList.get(i);

        this.missingFactorIndices = new ArrayList<>();
        for (int i = 0; i < numTaxa * dimTrait; ++i) {
            missingFactorIndices.add(i);
        }

        this.nuggetPrecision = nuggetPrecision;
        this.taskPool = (taskPool != null) ? taskPool : new TaskPool(numTaxa, 1);

        this.usePrecisionCache = cacheProvider.useCache();

        if (usePrecisionCache && this.taskPool.getNumThreads() > 1) {
            throw new IllegalArgumentException("Cannot currently parallelize cached precisions");
        }

        if (this.taskPool.getNumTaxon() != numTaxa) {
            throw new IllegalArgumentException("Incorrectly specified TaskPool");
        }

        this.observedInnerProduct = new DenseMatrix64F(numFactors, numFactors);


    }

    final private void setupObservedTraits(List<Integer> observedList, List<Integer> partialList) {
        for (int trait = 0; trait < dimTrait; trait++) {
            int nObserved = 0;
            for (int taxon = 0; taxon < numTaxa; taxon++) {
                if (observedIndicators[taxon][trait] == 1) {
                    nObserved += 1;
                }
            }
            if (nObserved == numTaxa) {
                observedList.add(trait);
            } else if (nObserved > 0) { // TODO: maybe change to `else {...` for path sampling
                partialList.add(trait);
            }
        }

    }

    @Override
    public boolean bufferTips() {
        return true;
    }

    @Override
    public int getTraitCount() {
        return 1;
    }

    @Override
    public int getDataDimension() {
        return dimTrait;
    }

    @Override
    public int getTraitDimension() {
        return numFactors;
    }  // Returns dimension of latent factors

    @Override
    public String getTipTraitName() {
        return tipTraitName;
    }

    @Override
    public void setTipTraitName(String name) {
        tipTraitName = name;
    }

    @Override
    public PrecisionType getPrecisionType() {
        return PrecisionType.FULL;
    }

    @Override
    public double[] getTipPartial(int taxonIndex, boolean fullyObserved) {
        if (fullyObserved) {
            throw new IllegalArgumentException("Wishart statistics are not implemented for the integrated factor model");
        }

        checkStatistics();

        double[] partial = new double[dimPartial];
        System.arraycopy(partials, taxonIndex * dimPartial, partial, 0, dimPartial);
        return partial;
    }

    @Override
    public List<Integer> getMissingIndices() {
        return missingFactorIndices;
    }

    @Override
    public boolean[] getDataMissingIndicators() {
        return missingDataIndicator;
    }

    @Override
    public boolean[] getTraitMissingIndicators() {
        if (getDataMissingIndicators() == null) {
            return null;
        } else if (missingTraitIndicators == null) {
            this.missingTraitIndicators = new boolean[getParameter().getDimension()];
            Arrays.fill(missingTraitIndicators, true); // all traits are latent
        }
        return missingTraitIndicators;
    }

    public List<Integer> getMissingDataIndices() {
        return missingDataIndices;
    }

    @Override
    public CompoundParameter getParameter() {
        return traitParameter;
    }

    @Override
    public boolean usesMissingIndices() {
        return true;
    }

    @Override
    public ContinuousTraitPartialsProvider[] getChildModels() {
        return new ContinuousTraitPartialsProvider[0]; // LFM is not currently extendible
    }

    @Override
    public boolean getDefaultAllowSingular() {
        return true;
    }

    @Override
    public boolean suppliesWishartStatistics() {
        return false;
    }

    @Override
    public Model getModel() {
        return this;
    }

    @Override
    public double getLogLikelihood() {
        if (!likelihoodKnown) {
            logLikelihood = calculateLogLikelihood();
            likelihoodKnown = true;
        }
        return 0;
    }

    @Override
    public void makeDirty() {
        likelihoodKnown = false;
        statisticsKnown = false;
        innerProductsKnown = false;
        observedInnerProductKnown = false;
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        // No model dependencies
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        observedInnerProductKnown = false;

        if (variable == loadingsTransposed) {
            statisticsKnown = false;
            likelihoodKnown = false;
            fireModelChanged(this);
        } else if (variable == traitParameter || variable == traitPrecision) {
            innerProductsKnown = false; // TODO: why does this not go to false when the loadings change???
            statisticsKnown = false;
            likelihoodKnown = false;
            fireModelChanged(this);
        } else {
            throw new RuntimeException("Unhandled parameter change type");
        }
    }

    @Override
    protected void storeState() {
        storedLogLikelihood = logLikelihood;
        storedLikelihoodKnow = likelihoodKnown;
        storedStatisticsKnown = statisticsKnown;

        System.arraycopy(partials, 0, storedPartials, 0, partials.length);
        System.arraycopy(normalizationConstants, 0,
                storedNormalizationConstants, 0, normalizationConstants.length);

        if (USE_INNER_PRODUCT_CACHE) {
            storedInnerProductsKnown = innerProductsKnown;

            System.arraycopy(traitInnerProducts, 0,
                    storedTraitInnerProducts, 0, traitInnerProducts.length);
        }
    }

    @Override
    protected void restoreState() {
        logLikelihood = storedLogLikelihood;
        likelihoodKnown = storedLikelihoodKnow;
        statisticsKnown = storedStatisticsKnown;

        double[] tmp1 = partials;
        partials = storedPartials;
        storedPartials = tmp1;

        double[] tmp2 = normalizationConstants;
        normalizationConstants = storedNormalizationConstants;
        storedNormalizationConstants = tmp2;

        if (USE_INNER_PRODUCT_CACHE) {
            innerProductsKnown = storedInnerProductsKnown;

            double[] tmp3 = traitInnerProducts;
            traitInnerProducts = storedTraitInnerProducts;
            storedTraitInnerProducts = tmp3;
        }

        observedInnerProductKnown = false; // TODO: proper store/restore
    }

    @Override
    protected void acceptState() {
        // Do nothing
    }

    // Access for FactorAnalysisOperatorAdaptor

    public int getNumberOfFactors() {
        return numFactors;
    }

    public int getNumberOfTaxa() {
        return numTaxa;
    }

    public int getNumberOfTraits() {
        return dimTrait;
    }

    public MatrixParameterInterface getLoadings() {
        return loadingsTransposed;
    }

    public Parameter getPrecision() {
        return traitPrecision;
    }

    // Private class functions

    private double calculateLogLikelihood() {

        checkStatistics();

        double logLikelihood = 0.0;
        for (double r : normalizationConstants) {
            logLikelihood += r;
        }
        return logLikelihood;
    }

    private void setupStatistics() {
        if (partials == null) {
            partials = new double[numTaxa * dimPartial];
            storedPartials = new double[numTaxa * dimPartial];
        }

        if (normalizationConstants == null) {
            normalizationConstants = new double[numTaxa];
            storedNormalizationConstants = new double[numTaxa];
        }

        if (USE_INNER_PRODUCT_CACHE) {
            if (traitInnerProducts == null) {
                traitInnerProducts = new double[numTaxa];
                storedTraitInnerProducts = new double[numTaxa];
            }
        }

        if (USE_INNER_PRODUCT_CACHE) {
            if (!innerProductsKnown) {
                setupInnerProducts();
                innerProductsKnown = true;
            }
        }

        loadings = loadingsTransposed.getParameterValues();
        gamma = traitPrecision.getParameterValues();

        computePartialsAndRemainders();
    }

    private final double nuggetPrecision;
    private final TaskPool taskPool;

    @Override
    public ContinuousExtensionDelegate getExtensionDelegate(ContinuousDataLikelihoodDelegate delegate,
                                                            TreeTrait treeTrait, Tree tree) {
        return new ContinuousExtensionDelegate.IndependentNormalExtensionDelegate(delegate, treeTrait,
                this, tree);
    }

    @Override
    public boolean diagonalVariance() {
        return true;
    }

    @Override
    public DenseMatrix64F getExtensionVariance() { //TODO: setup buffer if needed (probably not)
        //TODO: check that this does what it's supposed to.
        double[] precisionBuffer = traitPrecision.getParameterValues();
        DenseMatrix64F varianceMat = MissingOps.wrapDiagonalInverse(precisionBuffer, 0, precisionBuffer.length);
        return varianceMat;
    }

    @Override
    public DenseMatrix64F getExtensionVariance(NodeRef node) {
        return getExtensionVariance();
    }

    @Override
    public MatrixParameterInterface getExtensionPrecision() {
        //TODO: check that this does what it's supposed to.
        return new DiagonalMatrix(traitPrecision);
    }

    @Override
    public double[] transformTreeTraits(double[] treeTraits) {
        //TODO: check that this does what it's supposed to.

        DenseMatrix64F treeTraitMatrix = DenseMatrix64F.wrap(numTaxa, numFactors, treeTraits);
        DenseMatrix64F loadingsMatrix = DenseMatrix64F.wrap(numFactors, dimTrait,
                loadingsTransposed.getParameterValues());


        DenseMatrix64F traitMatrix = new DenseMatrix64F(numTaxa, dimTrait);
        org.ejml.ops.CommonOps.mult(treeTraitMatrix, loadingsMatrix, traitMatrix);

        if (DEBUG) {
            treeTraitMatrix.print();
            loadingsMatrix.print();
            traitMatrix.print();
        }

        return traitMatrix.data;
    }

    @Override
    public void chainRuleWrtVariance(double[] gradient, NodeRef node) {
        throw new RuntimeException("not yet implemented");
    }

    private void computePrecisionForTaxon(final DenseMatrix64F precision, final int taxon,
                                          final int numFactors) {

        final double[] observed = observedIndicators[taxon]; // TODO: only store for partiallyMissing?

        HashedMissingArray observedArray = null;
        DenseMatrix64F hashedPrecision = null;

        if (usePrecisionCache) {
            observedArray = new HashedMissingArray(observed);
            hashedPrecision = precisionMatrixMap.get(observedArray);
        }


        if (!usePrecisionCache || hashedPrecision == null) { // TODO: remove code duplication with below
            if (!observedInnerProductKnown) {
                for (int row = 0; row < numFactors; ++row) {
                    for (int col = row; col < numFactors; ++col) {
                        double sum = 0;
                        for (int k : fullyObservedTraits) {

                            sum += loadings[row * dimTrait + k] * //loadingsTransposed.getParameterValue(k, row) *
                                    gamma[k] *
                                    loadings[col * dimTrait + k]; // loadingsTransposed.getParameterValue(k, col);
                        }
                        observedInnerProduct.set(row, col, sum);
                        observedInnerProduct.set(col, row, sum);

                    }
                }

                observedInnerProductKnown = true;
            }

            // Compute L D_i \Gamma D_i^t L^t
            for (int row = 0; row < numFactors; ++row) {
                for (int col = row; col < numFactors; ++col) {
                    double sum = observedInnerProduct.get(row, col);
                    for (int k : partiallyMissingTraits) {
                        double thisPrecision = (observed[k] == 1.0) ?
                                gamma[k] // traitPrecision.getParameterValue(k)
                                : nuggetPrecision;
                        sum += loadings[row * dimTrait + k] * //loadingsTransposed.getParameterValue(k, row) *
                                thisPrecision *
                                loadings[col * dimTrait + k]; // loadingsTransposed.getParameterValue(k, col);
                    }
                    precision.unsafe_set(row, col, sum);
                    precision.unsafe_set(col, row, sum); // Symmetric matrix
                }
            }

            if (usePrecisionCache) {
                precisionMatrixMap.put(observedArray, precision);
            }

        } else {
            System.arraycopy(hashedPrecision.getData(), 0,
                    precision.getData(), 0, numFactors * numFactors);
        }
    }

    private static final boolean TIMING = false;
    private static final boolean USE_INNER_PRODUCT_CACHE = true;

    private Map<HashedMissingArray, DenseMatrix64F> precisionMatrixMap = new HashMap<>();

    private void fillInMeanForTaxon(final WrappedVector output, final DenseMatrix64F precision,
                                    final int taxon) {

        final double[] observed = observedIndicators[taxon];
//        final Parameter Y = traitParameter.getParameter(taxon);

        // Solve for a value \mu_i s.t. P_i \mu_i = (L D_i Y_i)

        final double[] tmp = new double[numFactors];
        final double[] tmp2 = new double[numFactors];

        for (int factor = 0; factor < numFactors; ++factor) {
            double sum = 0;
            for (int k = 0; k < dimTrait; ++k) {
                sum += loadings[factor * dimTrait + k] * //loadingsTransposed.getParameterValue(k, factor) *  // TODO Maybe a memory access issue here?
                        observed[k] * gamma[k] * // traitPrecision.getParameterValue(k) *
                        data[taxon * dimTrait + k];
//                        Y.getParameterValue(k);
            }
            tmp[factor] = sum;
        }

        DenseMatrix64F B = DenseMatrix64F.wrap(numFactors, 1, tmp);
        DenseMatrix64F X = DenseMatrix64F.wrap(numFactors, 1, tmp2);

        safeSolve(precision, B, X, false);

        for (int row = 0; row < numFactors; ++row) {
            output.set(row, X.unsafe_get(row, 0));
        }

//        return ci;
    }

    private double computeTraitInnerProduct(final int taxon) {
        final double[] observed = observedIndicators[taxon];
        final Parameter Y = traitParameter.getParameter(taxon);

        // Compute Y_i^t D_i^t \Gamma D_i Y_i
        double sum = 0;
        for (int k = 0; k < dimTrait; ++k) {
            sum += Y.getParameterValue(k) * Y.getParameterValue(k) *
                    observed[k] * traitPrecision.getParameterValue(k);
        }
        return sum;
    }

    private void cacheTraitInnerProducts(final int taxon) {
        traitInnerProducts[taxon] = computeTraitInnerProduct(taxon);
    }

    private void setupInnerProducts() {

        data = traitParameter.getParameterValues();

        if (TIMING) {
            for (int taxon = 0; taxon < numTaxa; ++taxon) {
                cacheTraitInnerProducts(taxon);
            }
        } else {
            taskPool.fork((taxon, thread) -> cacheTraitInnerProducts(taxon));
        }
    }

    private double computeFactorInnerProduct(final WrappedVector mean, final DenseMatrix64F precision) {
        // Compute \mu_i^t P_i \mu^t
        double sum = 0;
        for (int row = 0; row < numFactors; ++row) {
            for (int col = 0; col < numFactors; ++col) {
                sum += mean.get(row) * precision.unsafe_get(row, col) * mean.get(col);
            }
        }
        return sum;
    }

    private double getTraitLogDeterminant(final int taxon) {

        final double[] observed = observedIndicators[taxon];

        // Compute det( D_i \Gamma D_i^t)
        double logDet = 0.0;
        for (int k = 0; k < dimTrait; ++k) {
            if (observed[k] == 1.0) {
                logDet += Math.log(traitPrecision.getParameterValue(k));
            }
        }
        return logDet;
    }

    private void makeCompletedUnobserved(final DenseMatrix64F matrix, double diagonal) {
        for (int row = 0; row < numFactors; ++row) {
            for (int col = 0; col < numFactors; ++col) {
                double x = (row == col) ? diagonal : 0.0;
                matrix.unsafe_set(row, col, x);
            }
        }
    }

    private void computePartialAndRemainderForOneTaxon(int taxon,
                                                       DenseMatrix64F precision,
                                                       DenseMatrix64F variance) {

        final int partialsOffset = dimPartial * taxon;
        // Work with mean in-place
        final WrappedVector mean = new WrappedVector.Raw(partials, partialsOffset, numFactors);

        computePrecisionForTaxon(precision, taxon, numFactors);
        fillInMeanForTaxon(mean, precision, taxon);

        if (DEBUG) {
            System.err.println("taxon " + taxon);
            System.err.println("\tprecision: " + precision);
        }

        double constant;
        double nuggetDensity = 0;

        int effDim = 0;
        double factorLogDeterminant = precisionType.getMissingDeterminantValue();

        if (observedDimensions[taxon] == 0) {

            makeCompletedUnobserved(precision, 0);
            makeCompletedUnobserved(variance, Double.POSITIVE_INFINITY);
            constant = 0.0;

        } else {


            if (DEBUG) {
                System.err.println("\tmean: " + mean);
                //System.err.println("\n");
            }


            InversionResult ci = safeDeterminant(precision, false); //TODO: figure out how to remove this (I don't want to do it twice) (see safeMultivariateIntegrator.IncreaseVariances)
            effDim = ci.getEffectiveDimension();
            factorLogDeterminant = ci.getReturnCode() == InversionResult.Code.NOT_OBSERVED ? 0 : ci.getLogDeterminant();
//            factorLogDeterminant = ci.getLogDeterminant();
            double traitLogDeterminant = getTraitLogDeterminant(taxon);

//            final double logDetChange = traitLogDeterminant - factorLogDeterminant;

            final double factorInnerProduct = computeFactorInnerProduct(mean, precision);
            final double traitInnerProduct = USE_INNER_PRODUCT_CACHE ?
                    traitInnerProducts[taxon] : computeTraitInnerProduct(taxon);
            final double innerProductChange = traitInnerProduct - factorInnerProduct;

//            int dimensionChange = observedDimensions[taxon] - ci.getEffectiveDimension(); //TODO: use this effective dimension in safeMultivariateIntegrator

            if (DEBUG) {
                System.err.println("fIP: " + factorInnerProduct);
                System.err.println("tIP: " + traitInnerProduct);
//                System.err.println("fDet: " + factorLogDeterminant);
                System.err.println("tDet: " + traitLogDeterminant);
//                System.err.println("deltaDim: " + dimensionChange)
                System.err.println(" deltaIP: " + innerProductChange + "\n\n");
            }

//            constant = 0.5 * (logDetChange - innerProductChange) - LOG_SQRT_2_PI * (dimensionChange) -
//                    nuggetDensity;

            constant = 0.5 * (traitLogDeterminant - factorLogDeterminant - innerProductChange) -
                    LOG_SQRT_2_PI * (observedDimensions[taxon] - effDim) -
                    nuggetDensity;

        }

        // store in precision, variance and normalization constant
        unwrap(precision, partials, partialsOffset + numFactors); //TODO: use PrecisionType.fillPrecisionInPartials()
        precisionType.fillEffDimInPartials(partials, partialsOffset, effDim, numFactors);
        precisionType.fillDeterminantInPartials(partials, partialsOffset, factorLogDeterminant, numFactors);
        precisionType.fillRemainderInPartials(partials, partialsOffset, constant, numFactors);

        if (STORE_VARIANCE) {
            safeInvert2(precision, variance, true);
            unwrap(variance, partials, partialsOffset + numFactors + numFactors * numFactors);
        }

        normalizationConstants[taxon] = constant;
    }

    private void computePartialsAndRemainders() {

        final DenseMatrix64F[] precisions = new DenseMatrix64F[taskPool.getNumThreads()];
        final DenseMatrix64F[] variances = new DenseMatrix64F[taskPool.getNumThreads()];

        for (int i = 0; i < taskPool.getNumThreads(); ++i) {
            precisions[i] = new DenseMatrix64F(numFactors, numFactors);
            variances[i] = new DenseMatrix64F(numFactors, numFactors);
        }

        if (usePrecisionCache) {
            precisionMatrixMap.clear();
            if (DEBUG) {
                System.err.println("Hash CLEARED");
            }
        }

        if (TIMING) { // Do not use threads or lambda when timing
            for (int taxon = 0; taxon < numTaxa; ++taxon) {
                computePartialAndRemainderForOneTaxon(taxon, precisions[0], variances[0]);
            }
        } else {
            taskPool.fork((taxon, thread) ->
                    computePartialAndRemainderForOneTaxon(taxon, precisions[thread], variances[thread]));
        }
    }

    private static final boolean STORE_VARIANCE = true;
    private static final boolean DEBUG = false;

    private void checkStatistics() {
        synchronized (this) {
            if (!statisticsKnown) {
                setupStatistics();
                statisticsKnown = true;
            }
        }
    }

    private static double[][] setupObservedIndicators(List<Integer> missingIndices, int nTaxa, int dimTrait) {
        double[][] observed = new double[nTaxa][dimTrait];

        for (double[] v : observed) {
            Arrays.fill(v, 1.0);
        }

        if (missingIndices != null) {
            for (Integer idx : missingIndices) {
                int taxon = idx / dimTrait;
                int trait = idx % dimTrait;
                observed[taxon][trait] = 0.0;
            }
        }

        return observed;
    }

    private static int[] setupObservedDimensions(double[][] observed) {
        int length = observed.length;
        int[] dimensions = new int[length];

        for (int i = 0; i < length; ++i) {
            double sum = 0;
            for (double x : observed[i]) {
                sum += x;
            }
            dimensions[i] = (int) sum;
        }

        return dimensions;
    }

    private boolean likelihoodKnown = false;
    private boolean storedLikelihoodKnow;

    private boolean statisticsKnown = false;
    private boolean storedStatisticsKnown;

    private boolean innerProductsKnown = false;
    private boolean storedInnerProductsKnown;

    private double logLikelihood;
    private double storedLogLikelihood;

    private double[] partials;
    private double[] storedPartials;

    private double[] normalizationConstants;
    private double[] storedNormalizationConstants;

    private double[] traitInnerProducts;
    private double[] storedTraitInnerProducts;

    private double[] data;
    private double[] loadings;
    private double[] gamma;

    private final int numTaxa;
    private final int dimTrait;
    private final int dimPartial;
    private final int numFactors;
    private final CompoundParameter traitParameter;
    private final MatrixParameterInterface loadingsTransposed;
    private final Parameter traitPrecision;
    private final List<Integer> missingFactorIndices;
    private final List<Integer> missingDataIndices;
    private final boolean[] missingDataIndicator;

    private final double[][] observedIndicators;
    private final int[] observedDimensions;

    private final boolean usePrecisionCache;


    private static double LOG_SQRT_2_PI = 0.5 * Math.log(2 * Math.PI);

    //TODO: remove code duplicaton?
    public enum CacheProvider {
        USE_CACHE {
            @Override
            boolean useCache() {
                return true;
            }

        },
        NO_CACHE {
            @Override
            boolean useCache() {
                return false;
            }
        };

        abstract boolean useCache();

    }


    // TODO Move remainder into separate class file
    public static AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {
        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            MutableTreeModel treeModel = (MutableTreeModel) xo.getChild(MutableTreeModel.class);
            TreeTraitParserUtilities utilities = new TreeTraitParserUtilities();

            TreeTraitParserUtilities.TraitsAndMissingIndices returnValue =
                    utilities.parseTraitsFromTaxonAttributes(xo,
                            treeModel, true);
            CompoundParameter traitParameter = returnValue.traitParameter;
            boolean[] missingIndicators = returnValue.getMissingIndicators();

            MatrixParameterInterface loadings = (MatrixParameterInterface) xo.getElementFirstChild(LOADINGS);
            Parameter traitPrecision = (Parameter) xo.getElementFirstChild(PRECISION);

            double nugget = xo.getAttribute(NUGGET, 0.0);

            TaskPool taskPool = (TaskPool) xo.getChild(TaskPool.class);

            CacheProvider cacheProvider;
            boolean useCache = xo.getAttribute(CACHE_PRECISION, false);
            if (useCache) {
                cacheProvider = CacheProvider.USE_CACHE;
            } else {
                cacheProvider = CacheProvider.NO_CACHE;
            }


            return new IntegratedFactorAnalysisLikelihood(xo.getId(), traitParameter, missingIndicators,
                    loadings, traitPrecision, nugget, taskPool, cacheProvider);
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        @Override
        public String getParserDescription() {
            return null;
        }

        @Override
        public Class getReturnType() {
            return IntegratedFactorAnalysisLikelihood.class;
        }

        @Override
        public String getParserName() {
            return INTEGRATED_FACTOR_Model;
        }
    };

    public static final String INTEGRATED_FACTOR_Model = "integratedFactorModel";
    private static final String LOADINGS = "loadings";
    private static final String PRECISION = "precision";
    private static final String NUGGET = "nugget";
    private static final String CACHE_PRECISION = "cachePrecision";

    private final static XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(LOADINGS, new XMLSyntaxRule[]{
                    new ElementRule(MatrixParameterInterface.class),
            }),
            new ElementRule(PRECISION, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class),
            }),
            // Tree trait parser
            new ElementRule(MutableTreeModel.class),
            AttributeRule.newStringRule(TreeTraitParserUtilities.TRAIT_NAME),
            new ElementRule(TreeTraitParserUtilities.TRAIT_PARAMETER, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }),
            new ElementRule(TreeTraitParserUtilities.MISSING, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }, true),
            AttributeRule.newDoubleRule(NUGGET, true),
            AttributeRule.newBooleanRule(STANDARDIZE, true),
            new ElementRule(TaskPool.class, true),
            AttributeRule.newDoubleRule(TARGET_SD, true),
            AttributeRule.newBooleanRule(CACHE_PRECISION, true)

    };

    public void setLikelihoodDelegate(ContinuousDataLikelihoodDelegate delegate) {
        this.delegate = delegate;
    }

    private ContinuousDataLikelihoodDelegate delegate = null;

    private static Matrix buildDiagonalMatrix(double[] diagonals) {
        Matrix mat = new Matrix(diagonals.length, diagonals.length);
        for (int i = 0; i < diagonals.length; ++i) {
            mat.set(i, i, diagonals[i]);
        }
        return mat;
    }

    @Override
    public String getReport() {

        StringBuilder sb = new StringBuilder();

        double logComponents = 0;

        if (delegate != null) {

            double logInc = delegate.getCallbackLikelihood().getLogLikelihood();

            final Tree tree = delegate.getCallbackLikelihood().getTree();
            final BranchRates branchRates = delegate.getCallbackLikelihood().getBranchRateModel();
            sb.append(tree.toString());
            sb.append("\n\n");

            final double normalization = delegate.getRateTransformation().getNormalization();
            final double priorSampleSize = delegate.getRootProcessDelegate().getPseudoObservations();

            double[][] treeStructure = MultivariateTraitDebugUtilities.getTreeVariance(tree, branchRates, 1.0, Double.POSITIVE_INFINITY);
            sb.append("Tree structure:\n");
            sb.append(new Matrix(treeStructure));
            sb.append("\n\n");

            double[][] treeSharedLengths = MultivariateTraitDebugUtilities.getTreeVariance(tree, branchRates, normalization, Double.POSITIVE_INFINITY);

            double[][] treeVariance = MultivariateTraitDebugUtilities.getTreeVariance(tree, branchRates, normalization, priorSampleSize);

            Matrix treeV = new Matrix(treeVariance);
            Matrix treeP = treeV.inverse();

            sb.append("Tree variance:\n");
            sb.append(treeV);
            sb.append("Tree precision:\n");
            sb.append(treeP);
            sb.append("\n\n");

            Matrix Lt = new Matrix(loadingsTransposed.getParameterAsMatrix());
            sb.append("Loadings:\n");
            sb.append(Lt.transpose());
            sb.append("\n\n");

            double[][] diffusionPrecision = delegate.getDiffusionModel().getPrecisionmatrix();
            Matrix diffusionVariance = new Matrix(diffusionPrecision).inverse();

            Matrix loadingsVariance = null;
            try {
                loadingsVariance = Lt.product(diffusionVariance.product(Lt.transpose()));
            } catch (IllegalDimension illegalDimension) {
                illegalDimension.printStackTrace();
            }
            sb.append("Loadings variance:\n");
            sb.append(loadingsVariance);
            sb.append("\n\n");

            assert (loadingsVariance != null);

            Matrix loadingsFactorsVariance = MultivariateTraitDebugUtilities.getJointVarianceFactor(priorSampleSize,
                    treeVariance, treeSharedLengths, loadingsVariance.toComponents(), diffusionVariance.toComponents(),
                    delegate.getDiffusionProcessDelegate(), Lt);

            Matrix gamma = buildDiagonalMatrix(traitPrecision.getParameterValues());
            sb.append("Trait precision:\n");
            sb.append(gamma);
            sb.append("\n\n");
            Matrix gammaVariance = gamma.inverse();

            double[] tmp = new double[tree.getExternalNodeCount()];
            Arrays.fill(tmp, 1.0);
            Matrix identity = buildDiagonalMatrix(tmp);
            Matrix errorVariance = new Matrix(KroneckerOperation.product(identity.toComponents(), gammaVariance.toComponents()));

            sb.append("Loadings-factors variance:\n");
            sb.append(loadingsFactorsVariance);
            sb.append("\n\n");

            sb.append("Error variance\n");
            sb.append(errorVariance);
            sb.append("\n\n");

            Matrix totalVariance = null;
            try {
                totalVariance = loadingsFactorsVariance.add(errorVariance);
            } catch (IllegalDimension illegalDimension) {
                illegalDimension.printStackTrace();
            }

            double[] allData = getParameter().getParameterValues();

            List<Integer> notMissing = new ArrayList<>();
            for (int taxon = 0; taxon < numTaxa; ++taxon) {
                double[] observed = observedIndicators[taxon];
                for (int trait = 0; trait < dimTrait; ++trait) {
                    if (observed[trait] == 0.0) {
                        System.err.println("Missing taxon " + taxon + " trait " + trait);
                    } else {
                        notMissing.add(taxon * dimTrait + trait);
                    }
                }
            }

            double[] priorMean = delegate.getRootPrior().getMean();
            Matrix treeDrift = new Matrix(MultivariateTraitDebugUtilities.getTreeDrift(tree, priorMean, delegate.getIntegrator(), delegate.getDiffusionProcessDelegate()));

            if (delegate.getDiffusionProcessDelegate().hasDrift()) {
                sb.append("Tree drift (including root mean):\n");
                sb.append(new Matrix(treeDrift.toComponents()));
                sb.append("\n\n");
            }

            try {
                loadingsFactorsVariance = treeDrift.product(Lt.transpose());
            } catch (IllegalDimension illegalDimension) {
                illegalDimension.printStackTrace();
            }
            double[] drift = KroneckerOperation.vectorize(loadingsFactorsVariance.toComponents());

            int[] notMissingIndices = new int[notMissing.size()];
            double[] data = new double[notMissing.size()];
            for (int i = 0; i < notMissing.size(); ++i) {
                notMissingIndices[i] = notMissing.get(i);
                data[i] = allData[notMissing.get(i)];
            }

            if (totalVariance != null) {
                totalVariance = new Matrix(Matrix.gatherRowsAndColumns(totalVariance.toComponents(), notMissingIndices, notMissingIndices));
            }
            Matrix totalPrecision = null;
            if (totalVariance != null) {
                totalPrecision = totalVariance.inverse();
            }

            drift = Matrix.gatherEntries(drift, notMissingIndices);

            sb.append("Total variance:\n");
            sb.append(totalVariance);
            sb.append("\n\n");
            sb.append("Total precision:\n");
            sb.append(totalPrecision);
            sb.append("\n\n");

            sb.append("Data:\n");
            sb.append(new Vector(data));
            sb.append("\n\n");

            sb.append("Expectations:\n");
            sb.append(new Vector(drift));
            sb.append("\n\n");

            MultivariateNormalDistribution mvn = null;
            if (totalPrecision != null) {
                mvn = new MultivariateNormalDistribution(drift,
                        totalPrecision.toComponents());
            }

            double logDensity = 0;
            if (mvn != null) {
                logDensity = mvn.logPdf(data);
            }
            sb.append("logMultiVariateNormalDensity = ").append(logDensity).append("\n\n");

            sb.append("traitDataLikelihood = ").append(logInc).append("\n");
            logComponents += logInc;
        }

        sb.append("logLikelihood = ").append(getLogLikelihood()).append("\n");

        if (logComponents != 0.0) {
            sb.append("total likelihood = ").append((getLogLikelihood() + logComponents)).append("\n");
        }

        return sb.toString();
    }
}
