/*
 * IntegratedFactorTraitDataModel.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.evolution.tree.MultivariateTraitTree;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.continuous.cdi.PrecisionType;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.inference.model.*;
import dr.math.matrixAlgebra.missingData.InversionResult;
import dr.xml.*;
import org.ejml.data.DenseMatrix64F;

import java.util.Arrays;
import java.util.List;

import static dr.math.matrixAlgebra.missingData.MissingOps.safeInvert;
import static dr.math.matrixAlgebra.missingData.MissingOps.unwrap;

/**
 * @author Marc A. Suchard
 */

public class IntegratedFactorAnalysisLikelihood extends AbstractModelLikelihood
        implements ContinuousTraitPartialsProvider {

    public IntegratedFactorAnalysisLikelihood(String name,
                                              CompoundParameter traitParameter,
                                              List<Integer> missingIndices,
                                              MatrixParameterInterface loadings,
                                              Parameter traitPrecision) {
        super(name);

        this.traitParameter = traitParameter;
        this.loadings = loadings;
        this.traitPrecision = traitPrecision; // TODO Generalize for non-diagonal precision
//        this.missingIndices = missingIndices;

        this.numTaxa = traitParameter.getParameterCount();
        this.dimTrait = traitParameter.getParameter(0).getDimension();
        this.dimPartial =  dimTrait + PrecisionType.FULL.getMatrixLength(dimTrait);
        this.numFactors = loadings.getRowDimension();

        assert(dimTrait == loadings.getColumnDimension());

        addVariable(traitParameter);
        addVariable(loadings);
        addVariable(traitPrecision);

        this.observedIndicators = setupObservedIndicators(missingIndices, numTaxa, dimTrait);
        this.observedDimensions = setupObservedDimensions(observedIndicators);
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
            throw new IllegalArgumentException("Unsure what these values should be"); // TODO
        }

        checkStatistics();

        double[] partial = new double[dimPartial];
        System.arraycopy(partials, taxonIndex * dimPartial, partial, 0, dimPartial);
        return partial;
    }

//    @Override
//    public double[] getTipPartial(int taxonIndex) {
//        checkStatistics();
//        throw new RuntimeException("To implement");
//    }

//    @Override
//    public double[] getTipObservation(int taxonIndex, PrecisionType precisionType) {
//        checkStatistics();
//        throw new RuntimeException("To implement");
//    }

    @Override
    public List<Integer> getMissingIndices() {
        return null; // TODO Fix use-case (all tree-tip values (factors) should be missing)
    }

    @Override
    public CompoundParameter getParameter() {
        return null; // TODO Fix use-case
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
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        // No model dependencies
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if (variable == loadings || variable == traitPrecision) {
            statisticsKnown = false;
            likelihoodKnown = false;
            fireModelChanged(this);
//            fireModelChanged(this, getTaxonIndex(index));
        } else if (variable == traitParameter) {
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
        System.arraycopy(remainders, 0, storedRemainders, 0, remainders.length);
    }

    @Override
    protected void restoreState() {
        logLikelihood = storedLogLikelihood;
        likelihoodKnown = storedLikelihoodKnow;
        statisticsKnown = storedStatisticsKnown;

        double[] tmp1 = partials;
        partials = storedPartials;
        storedPartials = tmp1;

        double[] tmp2 = remainders;
        remainders = storedRemainders;
        storedRemainders = tmp2;
    }

    @Override
    protected void acceptState() {
        // Do nothing
    }

    private double calculateLogLikelihood() {
        double logLikelihood = 0.0;
        for (double r : remainders) {
            logLikelihood += r;
        }
        return logLikelihood;
    }

    private void setupStatistics() {
        if (partials == null) {
            partials = new double[numTaxa * dimPartial];
            storedPartials = new double[numTaxa * dimPartial];
        }

        if (remainders == null) {
            remainders = new double[numTaxa];
            storedRemainders = new double[numTaxa];
        }

        computePartialsAndRemainders();
    }


    private void fillPrecisionForTaxon(final DenseMatrix64F precision, final int taxon) {

        final double[] observed = observedIndicators[taxon];

        // Compute L D_i \Gamma L^t   // TODO Generalize for non-diagonal \Gamma
        for (int row = 0; row < numFactors; ++row) {
            for (int col = 0; col < numFactors; ++col) {
                double sum = 0;
                for (int k = 0; k < dimTrait; ++k) {
                    sum += loadings.getParameterValue(row, k) *
                            observed[k] * traitPrecision.getParameterValue(k) *
                            loadings.getParameterValue(col, k);
                }
                precision.unsafe_set(row, col, sum);
            }
        }
    }

    private double getTraitDeterminant(final int taxon) {

        final double[] observed = observedIndicators[taxon];

        // Compute det( D_i \Gamma ) // TODO Generalize for non-diagonal \Gamma
        double det = 1.0;
        for (int k = 0; k < dimTrait; ++k) {
            if (observed[k] == 1.0) {
                det *= traitPrecision.getParameterValue(k);
            }
        }
        return det;
    }

    private void computePartialsAndRemainders() {
        
        final DenseMatrix64F precision = new DenseMatrix64F(dimTrait, dimTrait);
        final DenseMatrix64F variance = new DenseMatrix64F(dimTrait, dimTrait);

        int partialsOffset = 0;
        for (int taxon = 0; taxon < numTaxa; ++taxon) {

            fillPrecisionForTaxon(precision, taxon);
            InversionResult ci = safeInvert(precision, variance, true);

            final double factorDeterminant = ci.getDeterminant();
            final double traitDeterminant = getTraitDeterminant(taxon);

            // TODO Fill in partial mean starting at partials[partialsOffset]

            unwrap(precision, partials, partialsOffset + dimTrait);
            unwrap(variance, partials, partialsOffset + dimTrait + dimTrait * dimTrait);

            // TODO Fill in normalization constant into remainder[taxon]

            partialsOffset += dimPartial;
        }
    }

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

        for (Integer idx : missingIndices) {
            int taxon = idx / dimTrait;
            int trait = idx % dimTrait;
            observed[taxon][trait] = 0.0;
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
    private double logLikelihood;
    private double storedLogLikelihood;

    private double[] partials;
    private double[] storedPartials;

    private double[] remainders;
    private double[] storedRemainders;

    private final int numTaxa;
    private final int dimTrait;
    private final int dimPartial;
    private final int numFactors;
    private final CompoundParameter traitParameter;
    private final MatrixParameterInterface loadings;
    private final Parameter traitPrecision;
//    private final List<Integer> missingIndices;

    private final double[][] observedIndicators;
    private final int[] observedDimensions;

    // TODO Move remainder into separate class file
    public static AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {
        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            MultivariateTraitTree treeModel = (MultivariateTraitTree) xo.getChild(TreeModel.class);
            TreeTraitParserUtilities utilities = new TreeTraitParserUtilities();

            TreeTraitParserUtilities.TraitsAndMissingIndices returnValue =
                    utilities.parseTraitsFromTaxonAttributes(xo, TreeTraitParserUtilities.DEFAULT_TRAIT_NAME,
                            treeModel, true);
            CompoundParameter traitParameter = returnValue.traitParameter;
            List<Integer> missingIndices = returnValue.missingIndices;

            MatrixParameterInterface loadings = (MatrixParameterInterface) xo.getElementFirstChild(LOADINGS);
            Parameter traitPrecision = (Parameter) xo.getElementFirstChild(PRECISION);

            return new IntegratedFactorAnalysisLikelihood(xo.getId(), traitParameter, missingIndices,
                    loadings, traitPrecision);
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
    public static final String LOADINGS = "loadings";
    public static final String PRECISION = "precision";

    private final static XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
            new ElementRule(LOADINGS, new XMLSyntaxRule[] {
                    new ElementRule(MatrixParameterInterface.class),
            }),
            new ElementRule(PRECISION, new XMLSyntaxRule[] {
                    new ElementRule(Parameter.class),
            }),
            // Tree trait parser
            new ElementRule(MultivariateTraitTree.class),
            AttributeRule.newStringRule(TreeTraitParserUtilities.TRAIT_NAME),
            new ElementRule(TreeTraitParserUtilities.TRAIT_PARAMETER, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }),
            new ElementRule(TreeTraitParserUtilities.MISSING, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }, true),

    };
}
