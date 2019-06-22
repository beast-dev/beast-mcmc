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

import dr.evolution.tree.BranchRates;
import dr.evolution.tree.MutableTreeModel;
import dr.evolution.tree.Tree;

import dr.evomodel.continuous.hmc.TaxonTaskPool;
import dr.evomodel.treedatalikelihood.continuous.cdi.PrecisionType;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.inference.model.*;
import dr.math.KroneckerOperation;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.math.matrixAlgebra.IllegalDimension;
import dr.math.matrixAlgebra.Matrix;
import dr.math.matrixAlgebra.Vector;
import dr.math.matrixAlgebra.WrappedVector;
import dr.math.matrixAlgebra.missingData.InversionResult;
import dr.xml.*;
import org.ejml.data.DenseMatrix64F;

import java.util.*;

import static dr.evomodelxml.treelikelihood.TreeTraitParserUtilities.STANDARDIZE;
import static dr.evomodelxml.treelikelihood.TreeTraitParserUtilities.TARGET_SD;
import static dr.math.matrixAlgebra.missingData.MissingOps.safeInvert2;
import static dr.math.matrixAlgebra.missingData.MissingOps.safeSolve;
import static dr.math.matrixAlgebra.missingData.MissingOps.unwrap;

/**
 * @author Marc A. Suchard
 * @author Gabriel Hassler
 */

public class IntegratedFactorAnalysisLikelihood extends AbstractModelLikelihood
        implements ContinuousTraitPartialsProvider, Reportable {

    public IntegratedFactorAnalysisLikelihood(String name,
                                              CompoundParameter traitParameter,
                                              List<Integer> missingIndices,
                                              MatrixParameterInterface loadings,
                                              Parameter traitPrecision,
                                              double nuggetPrecision,
                                              TaxonTaskPool taxonTaskPool) {
        super(name);

        this.traitParameter = traitParameter;
        this.loadingsTransposed = loadings;
        this.traitPrecision = traitPrecision;

        this.numTaxa = traitParameter.getParameterCount();
        this.dimTrait = traitParameter.getParameter(0).getDimension();
        this.numFactors = loadings.getColumnDimension();

        assert (dimTrait == loadings.getRowDimension());

        this.dimPartial = numFactors + PrecisionType.FULL.getMatrixLength(numFactors);

        addVariable(traitParameter);
        addVariable(loadings);
        addVariable(traitPrecision);

        this.missingDataIndices = missingIndices;
        this.missingDataIndicator = ContinuousTraitPartialsProvider.indicesToIndicator(
                missingIndices, traitParameter.getDimension());
        this.observedIndicators = setupObservedIndicators(missingDataIndices, numTaxa, dimTrait);
        this.observedDimensions = setupObservedDimensions(observedIndicators);

        this.missingFactorIndices = new ArrayList<>();
        for (int i = 0; i < numTaxa * dimTrait; ++i) {
            missingFactorIndices.add(i);
        }

        this.nuggetPrecision = nuggetPrecision;
        this.taxonTaskPool = (taxonTaskPool != null) ? taxonTaskPool : new TaxonTaskPool(numTaxa, 1);

        if (USE_PRECISION_CACHE && this.taxonTaskPool.getNumThreads() > 1) {
            throw new IllegalArgumentException("Cannot currently parallelize cached precisions");
        }

        if (this.taxonTaskPool.getNumTaxon() != numTaxa) {
            throw new IllegalArgumentException("Incorrectly specified TaxonTaskPool");
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

    public int getDataDimension() {
        return dimTrait;
    }

    @Override
    public int getTraitDimension() {
        return numFactors;
    }  // Returns dimension of latent factors

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
    public boolean[] getMissingIndicator() {
        return missingDataIndicator;
    }

    public List<Integer> getMissingDataIndices() {
        return missingDataIndices;
    }

    @Override
    public CompoundParameter getParameter() {
        return traitParameter;
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
        return logLikelihood;
    }

    @Override
    public void makeDirty() {
        likelihoodKnown = false;
        statisticsKnown = false;
        innerProductsKnown = false;
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        // No model dependencies
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if (variable == loadingsTransposed) {
            statisticsKnown = false;
            likelihoodKnown = false;
            fireModelChanged(this);
        } else if (variable == traitParameter || variable == traitPrecision) {
            innerProductsKnown = false;
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

        computePartialsAndRemainders();
    }

    private final double nuggetPrecision;
    private final TaxonTaskPool taxonTaskPool;

    private class HashedMissingArray {

        final private double[] array;

        HashedMissingArray(final double[] array) {
            this.array = array;
        }

        public double[] getArray() {
            return array;
        }

        public double get(int index) {
            return array[index];
        }

        public int getLength() {
            return array.length;
        }

        @Override
        public int hashCode() {
            return Arrays.hashCode(array);
        }

        @Override
        public boolean equals(Object obj) {
            return obj instanceof HashedMissingArray && Arrays.equals(array,
                    ((HashedMissingArray) obj).array);
        }

        public String toString() {
            return new Vector(array).toString();
        }
    }

    private void computePrecisionForTaxon(final DenseMatrix64F precision, final int taxon,
                                          final int numFactors) {

        final double[] observed = observedIndicators[taxon];

        final HashedMissingArray observedArray;
        DenseMatrix64F hashedPrecision;

        if (USE_PRECISION_CACHE) {
            observedArray = new HashedMissingArray(observed);
            hashedPrecision = precisionMatrixMap.get(observedArray);
        }

        // TODO Only need to compute for each unique set of observed[] << numTaxa

        if (!USE_PRECISION_CACHE || hashedPrecision == null) {

            // Compute L D_i \Gamma D_i^t L^t
            for (int row = 0; row < numFactors; ++row) {
                for (int col = 0; col < numFactors; ++col) {
                    double sum = 0;
                    for (int k = 0; k < dimTrait; ++k) {
                        double thisPrecision = (observed[k] == 1.0) ?
                                traitPrecision.getParameterValue(k) : nuggetPrecision;
                        sum += loadingsTransposed.getParameterValue(k, row) *
                                thisPrecision *
                                loadingsTransposed.getParameterValue(k, col);
                    }
                    precision.unsafe_set(row, col, sum);
                }
            }

            if (USE_PRECISION_CACHE) {
                precisionMatrixMap.put(observedArray, precision);
            }

        } else {
            System.arraycopy(hashedPrecision.getData(), 0,
                    precision.getData(), 0, numFactors * numFactors);
        }
    }

    private static final boolean TIMING = false;
    private static final boolean USE_INNER_PRODUCT_CACHE = true;
    private static final boolean USE_PRECISION_CACHE = false;

    private Map<HashedMissingArray, DenseMatrix64F> precisionMatrixMap = new HashMap<>();

    private InversionResult fillInMeanForTaxon(final WrappedVector output, final DenseMatrix64F precision,
                                               final int taxon) {

        final double[] observed = observedIndicators[taxon];
        final Parameter Y = traitParameter.getParameter(taxon);

        // Solve for a value \mu_i s.t. P_i \mu_i = (L D_i Y_i)

        final double[] tmp = new double[numFactors];
        final double[] tmp2 = new double[numFactors];

        for (int row = 0; row < numFactors; ++row) {
            double sum = 0;
            for (int k = 0; k < dimTrait; ++k) {
                sum += loadingsTransposed.getParameterValue(k, row) *  // TODO Maybe a memory access issue here?
                        observed[k] * traitPrecision.getParameterValue(k) *
                        Y.getParameterValue(k);
            }
            tmp[row] = sum;
        }

        DenseMatrix64F B = DenseMatrix64F.wrap(numFactors, 1, tmp);
        DenseMatrix64F X = DenseMatrix64F.wrap(numFactors, 1, tmp2);

        InversionResult ci = safeSolve(precision, B, X, true);

        for (int row = 0; row < numFactors; ++row) {
            output.set(row, X.unsafe_get(row, 0));
        }

        return ci;
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
        if (TIMING) {
            for (int taxon = 0; taxon < numTaxa; ++taxon) {
                cacheTraitInnerProducts(taxon);
            }
        } else {
            taxonTaskPool.fork((taxon, thread) -> cacheTraitInnerProducts(taxon));
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
        InversionResult ci = fillInMeanForTaxon(mean, precision, taxon);

        if (DEBUG) {
            System.err.println("taxon " + taxon);
            System.err.println("\tprecision: " + precision);
        }

        double constant;
        double nuggetDensity = 0;

        if (observedDimensions[taxon] == 0) {

            makeCompletedUnobserved(precision, 0);
            makeCompletedUnobserved(variance, Double.POSITIVE_INFINITY);
            constant = 0.0;

        } else {


            if (DEBUG) {
                System.err.println("\tmean: " + mean);
                //System.err.println("\n");
            }

            final double factorLogDeterminant = ci.getLogDeterminant();
            double traitLogDeterminant = getTraitLogDeterminant(taxon);

            final double logDetChange = traitLogDeterminant - factorLogDeterminant;

            final double factorInnerProduct = computeFactorInnerProduct(mean, precision);
            final double traitInnerProduct = USE_INNER_PRODUCT_CACHE ?
                                traitInnerProducts[taxon] : computeTraitInnerProduct(taxon);
            final double innerProductChange = traitInnerProduct - factorInnerProduct;

            int dimensionChange = observedDimensions[taxon] - ci.getEffectiveDimension(); //TODO: use this effective dimension in safeMultivariateIntegrator

            if (DEBUG) {
                System.err.println("fIP: " + factorInnerProduct);
                System.err.println("tIP: " + traitInnerProduct);
                System.err.println("fDet: " + factorLogDeterminant);
                System.err.println("tDet: " + traitLogDeterminant);
                System.err.println("deltaDim: " + dimensionChange + " deltaIP: " + innerProductChange + "\n\n");
            }

            constant = 0.5 * (logDetChange - innerProductChange) - LOG_SQRT_2_PI * (dimensionChange) -
                    nuggetDensity;

        }

        // store in precision, variance and normalization constant
        unwrap(precision, partials, partialsOffset + numFactors);
        PrecisionType.FULL.fillEffDimInPartials(partials, partialsOffset, ci.getEffectiveDimension(), numFactors);

        if (STORE_VARIANCE) {
            safeInvert2(precision, variance, true);
            unwrap(variance, partials, partialsOffset + numFactors + numFactors * numFactors);
        }

        normalizationConstants[taxon] = constant;
    }

    private void computePartialsAndRemainders() {

        final DenseMatrix64F[] precisions = new DenseMatrix64F[taxonTaskPool.getNumThreads()];
        final DenseMatrix64F[] variances = new DenseMatrix64F[taxonTaskPool.getNumThreads()];

        for (int i = 0; i < taxonTaskPool.getNumThreads(); ++i) {
            precisions[i] = new DenseMatrix64F(numFactors, numFactors);
            variances[i] = new DenseMatrix64F(numFactors, numFactors);
        }

        if (USE_PRECISION_CACHE) {
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
            taxonTaskPool.fork((taxon, thread) ->
                    computePartialAndRemainderForOneTaxon(taxon, precisions[thread], variances[thread]));
        }
    }

    private static final boolean STORE_VARIANCE = true;
    private static final boolean DEBUG = false;

    private void checkStatistics() {
        if (!statisticsKnown) {
            setupStatistics();
            statisticsKnown = true;
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

    private static double LOG_SQRT_2_PI = 0.5 * Math.log(2 * Math.PI);

    // TODO Move remainder into separate class file
    public static AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {
        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            MutableTreeModel treeModel = (MutableTreeModel) xo.getChild(MutableTreeModel.class);
            TreeTraitParserUtilities utilities = new TreeTraitParserUtilities();

            TreeTraitParserUtilities.TraitsAndMissingIndices returnValue =
                    utilities.parseTraitsFromTaxonAttributes(xo, TreeTraitParserUtilities.DEFAULT_TRAIT_NAME,
                            treeModel, true);
            CompoundParameter traitParameter = returnValue.traitParameter;
            List<Integer> missingIndices = returnValue.missingIndices;

            MatrixParameterInterface loadings = (MatrixParameterInterface) xo.getElementFirstChild(LOADINGS);
            Parameter traitPrecision = (Parameter) xo.getElementFirstChild(PRECISION);

            double nugget = xo.getAttribute(NUGGET, 0.0);

            TaxonTaskPool taxonTaskPool = (TaxonTaskPool) xo.getChild(TaxonTaskPool.class);

            return new IntegratedFactorAnalysisLikelihood(xo.getId(), traitParameter, missingIndices,
                    loadings, traitPrecision, nugget, taxonTaskPool);
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

    private static final String INTEGRATED_FACTOR_Model = "integratedFactorModel";
    private static final String LOADINGS = "loadings";
    private static final String PRECISION = "precision";
    private static final String NUGGET = "nugget";

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
            new ElementRule(TaxonTaskPool.class, true),
            AttributeRule.newDoubleRule(TARGET_SD, true),

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
