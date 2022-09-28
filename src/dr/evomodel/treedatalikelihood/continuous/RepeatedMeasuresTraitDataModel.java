/*
 * RepeatedMeasuresTraitDataModel.java
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

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.treedatalikelihood.continuous.cdi.PrecisionType;
import dr.evomodel.treedatalikelihood.preorder.ContinuousExtensionDelegate;
import dr.evomodel.treedatalikelihood.preorder.ModelExtensionProvider;
import dr.evomodelxml.continuous.ContinuousTraitDataModelParser;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.inference.model.*;
import dr.math.matrixAlgebra.*;
import dr.math.matrixAlgebra.missingData.InversionResult;
import dr.math.matrixAlgebra.missingData.MissingOps;
import dr.xml.*;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import java.util.ArrayList;
import java.util.Arrays;

import static dr.evomodelxml.continuous.ContinuousTraitDataModelParser.NUM_TRAITS;

/**
 * @author Marc A. Suchard
 * @author Gabriel Hassler
 */
public class RepeatedMeasuresTraitDataModel extends ContinuousTraitDataModel implements FullPrecisionContinuousTraitPartialsProvider,
        ModelExtensionProvider.NormalExtensionProvider {

    private final String traitName;
    private final MatrixParameterInterface samplingPrecisionParameter;
    private boolean diagonalOnly = false;
    private boolean variableChanged = true;
    private boolean varianceKnown = false;

    private Matrix samplingPrecision;
    private Matrix samplingVariance;
    private Matrix storedSamplingPrecision;
    private Matrix storedSamplingVariance;
    private boolean storedVarianceKnown = false;
    private boolean storedVariableChanged = true;

    private boolean[] missingTraitIndicators = null;

    private ContinuousTraitPartialsProvider childModel;

    private final int nRepeats;
    private ArrayList<Integer>[] relevantRepeats;
    private final int nObservedTips;

    private final static double LOG2PI = Math.log(Math.PI * 2);


    public RepeatedMeasuresTraitDataModel(String name,
                                          ContinuousTraitPartialsProvider childModel,
                                          CompoundParameter parameter,
                                          boolean[] missindIndicators,
                                          boolean useMissingIndices,
                                          final int dimTrait,
                                          final int numTraits,
                                          MatrixParameterInterface samplingPrecision,
                                          PrecisionType precisionType) {

        super(name, parameter, missindIndicators, useMissingIndices, dimTrait, numTraits, precisionType);
        if (numTraits > 1) {
            throw new RuntimeException("not currently implemented");
        }

        this.childModel = childModel;
        this.traitName = name;
        this.samplingPrecisionParameter = samplingPrecision;
        this.nRepeats = childModel.getTraitCount() / numTraits;
        addVariable(samplingPrecision);

        calculatePrecisionInfo();

//        this.samplingVariance = new Matrix(samplingPrecision.getParameterAsMatrix()).inverse();
        this.samplingVariance = null;


        samplingPrecisionParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0,
                samplingPrecision.getDimension()));


        int offsetInc = precisionType.getPartialsDimension(this.dimTrait);

        int nTaxa = getParameter().getParameterCount();
        int nObservedTips = 0;
        relevantRepeats = new ArrayList[nTaxa];

        for (int i = 0; i < nTaxa; i++) {
            int precisionOffset = precisionType.getPrecisionOffset(this.dimTrait);

            relevantRepeats[i] = new ArrayList<>();
            double[] partial = childModel.getTipPartial(i, false);
            for (int r = 0; r < nRepeats; r++) {
                boolean isAtLeastPartiallyObserved = false;
                DenseMatrix64F P = MissingOps.wrap(partial, precisionOffset, this.dimTrait, this.dimTrait);

                for (int j = 0; j < this.dimTrait; j++) {
                    if (P.get(j, j) > 0) {
                        isAtLeastPartiallyObserved = true;
                        break;
                    }
                }

                if (isAtLeastPartiallyObserved) {
                    relevantRepeats[i].add(r);
                    nObservedTips++;
                }

                precisionOffset += offsetInc;

            }
        }

        this.nObservedTips = nObservedTips;
    }

    @Override
    public double[] getTipPartial(int taxonIndex, boolean fullyObserved) {

        assert (numTraits == 1);
        assert (samplingPrecision.rows() == dimTrait && samplingPrecision.columns() == dimTrait);

        recomputeVariance();

        if (fullyObserved) {
            throw new RuntimeException("Incompatible with this model.");
        }

        double[] partial = childModel.getTipPartial(taxonIndex, fullyObserved);

        if (nRepeats == 1) {
            if (precisionType == precisionType.SCALAR) {
                return partial; //TODO: I don't think this is right, especially given constructor above.
            }
            DenseMatrix64F V = MissingOps.wrap(partial, dimTrait + dimTrait * dimTrait, dimTrait, dimTrait);

            //TODO: remove diagonalOnly part
            if (diagonalOnly) {
                for (int index = 0; index < dimTrait; index++) {
                    V.set(index, index, V.get(index, index) + 1 / samplingPrecision.component(index, index));
                }
            } else {
                for (int i = 0; i < dimTrait; i++) {
                    for (int j = 0; j < dimTrait; j++) {
                        V.set(i, j, V.get(i, j) + samplingVariance.component(i, j));
                    }
                }
            }


            DenseMatrix64F P = new DenseMatrix64F(dimTrait, dimTrait);
            MissingOps.safeInvert2(V, P, false); //TODO this isn't necessary when this is fully observed

            MissingOps.unwrap(P, partial, dimTrait);
            MissingOps.unwrap(V, partial, dimTrait + dimTrait * dimTrait);

            if (DEBUG) {
                System.err.println("taxon " + taxonIndex);
                System.err.println("\tprecision: " + P);
                System.err.println("\tmean: " + new WrappedVector.Raw(partial, 0, dimTrait));
            }

            return partial;
        }

        int offsetInc = precisionType.getPartialsDimension(dimTrait);
        int varOffset = precisionType.getVarianceOffset(dimTrait);
        int meanOffset = precisionType.getMeanOffset(dimTrait);
        int varDim = precisionType.getVarianceLength(dimTrait);
        int remOffset = precisionType.getRemainderOffset(dimTrait);

        DenseMatrix64F Pi = new DenseMatrix64F(dimTrait, dimTrait);
        DenseMatrix64F Vi = new DenseMatrix64F(dimTrait, dimTrait);
        DenseMatrix64F P = new DenseMatrix64F(dimTrait, dimTrait);
        DenseMatrix64F V = new DenseMatrix64F(dimTrait, dimTrait);
        DenseMatrix64F Pm = new DenseMatrix64F(dimTrait, 1);
        DenseMatrix64F m = new DenseMatrix64F(dimTrait, 1);

        double remainder = 0;


        for (int i : relevantRepeats[taxonIndex]) {

            System.arraycopy(partial, offsetInc * i + varOffset, Vi.data, 0, varDim);
            for (int row = 0; row < dimTrait; row++) {
                if (Vi.get(row, row) < Double.POSITIVE_INFINITY) {
                    Vi.set(row, row, Vi.get(row, row) + samplingVariance.component(row, row));
                    for (int col = 0; col < row; col++) {
                        if (Vi.get(col, col) < Double.POSITIVE_INFINITY) {
                            Vi.set(row, col, Vi.get(row, col) + samplingVariance.component(row, col));
                            Vi.set(col, row, Vi.get(row, col));
                        }
                    }
                }
            }

            InversionResult result = MissingOps.safeInvert2(Vi, Pi, true);

            CommonOps.addEquals(P, Pi);

//            System.arraycopy(partial, meanOffset, mi.data, 0, dimTrait);
            double sumSquares = 0;

            for (int row = 0; row < dimTrait; row++) {

                int offset = offsetInc * i + meanOffset;
                double mr = partial[offset + row];

                double value = 0;
                for (int col = 0; col < dimTrait; col++) {
                    double mc = partial[offset + col];
                    double x = Pi.get(row, col) * mc;
                    value += x;
                    sumSquares += x * mr;
                }
                Pm.add(row, 0, value);
            }

            remainder += partial[offsetInc * i + remOffset];

            remainder -= result.getEffectiveDimension() * LOG2PI + sumSquares + result.getLogDeterminant();

        }


        MissingOps.safeSolve(P, Pm, m, false);
        InversionResult result = MissingOps.safeInvertPrecision(P, V, true); //TODO: don't invert twice
        if (result.getReturnCode() == InversionResult.Code.NOT_OBSERVED) {
            remainder = 0;
        } else {
            double sumSquares = 0;
            for (int row = 0; row < dimTrait; row++) {
                for (int col = 0; col < dimTrait; col++) {
                    sumSquares += m.get(row, 0) * m.get(col, 0) * P.get(row, col);
                }
            }

            remainder += result.getEffectiveDimension() * LOG2PI + sumSquares - result.getLogDeterminant();
        }


        partial = new double[offsetInc];

        System.arraycopy(m.data, 0, partial, precisionType.getMeanOffset(dimTrait), dimTrait);
        System.arraycopy(P.data, 0, partial, precisionType.getPrecisionOffset(dimTrait), varDim);
        System.arraycopy(V.data, 0, partial, precisionType.getVarianceOffset(dimTrait), varDim);
        precisionType.fillRemainderInPartials(partial, 0, 0.5 * remainder, dimTrait);

        return partial;
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


    private void recomputeVariance() {
        checkVariableChanged();
        if (!varianceKnown) {
            samplingVariance = samplingPrecision.inverse();
            varianceKnown = true;
        }
    }

    public Matrix getSamplingVariance() {
        recomputeVariance();
        return samplingVariance;
    }

    public String getTraitName() {
        return traitName;
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        super.handleVariableChangedEvent(variable, index, type);

        if (variable == samplingPrecisionParameter) {

            variableChanged = true;
            varianceKnown = false;
            fireModelChanged();
        }
    }

    private void calculatePrecisionInfo() {
        samplingPrecision = new Matrix(samplingPrecisionParameter.getParameterAsMatrix());
    }

    private void checkVariableChanged() {
        if (variableChanged) {
            calculatePrecisionInfo();
            variableChanged = false;
            varianceKnown = false;
        }
    }

    @Override
    protected void storeState() {
        storedSamplingPrecision = samplingPrecision.clone();
        storedSamplingVariance = samplingVariance.clone();
        storedVarianceKnown = varianceKnown;
        storedVariableChanged = variableChanged;
    }

    @Override
    protected void restoreState() {
        Matrix tmp = samplingPrecision;
        samplingPrecision = storedSamplingPrecision;
        storedSamplingPrecision = tmp;

        tmp = samplingVariance;
        samplingVariance = storedSamplingVariance;
        storedSamplingVariance = tmp;

        varianceKnown = storedVarianceKnown;
        variableChanged = storedVariableChanged;
    }

    @Override
    public ContinuousExtensionDelegate getExtensionDelegate(ContinuousDataLikelihoodDelegate delegate,
                                                            TreeTrait treeTrait, Tree tree) {
        checkVariableChanged();
        return new ContinuousExtensionDelegate.MultivariateNormalExtensionDelegate(delegate, treeTrait,
                this, tree);
    }

    @Override
    public boolean diagonalVariance() {
        return false; //TODO: base on precisionType
    }

    @Override
    public DenseMatrix64F getExtensionVariance() {
        recomputeVariance();
        double[] buffer = samplingVariance.toArrayComponents();
        return DenseMatrix64F.wrap(dimTrait, dimTrait, buffer);
    }

    @Override
    public DenseMatrix64F getExtensionVariance(NodeRef node) {
        return getExtensionVariance();
    }

    @Override
    public MatrixParameterInterface getExtensionPrecision() {
        return getExtensionPrecisionParameter(); //TODO: deprecate
    }

    public void getMeanTipVariances(DenseMatrix64F samplingVariance, DenseMatrix64F samplingComponent) {
        CommonOps.scale(1.0, samplingVariance, samplingComponent);
    }

    @Override
    public MatrixParameterInterface getExtensionPrecisionParameter() {
        checkVariableChanged();
        return samplingPrecisionParameter;
    }

    @Override
    public int getDataDimension() {
        return dimTrait;
    }

    @Override
    public boolean suppliesWishartStatistics() {
        return false;
    }

    @Override
    public void chainRuleWrtVariance(double[] gradient, NodeRef node) {
        // Do nothing
    }

    @Override
    public ContinuousTraitPartialsProvider[] getChildModels() {
        return new ContinuousTraitPartialsProvider[]{childModel};
    }

    @Override
    public double[] drawTraitsBelowConditionalOnDataAndTraitsAbove(double[] aboveTraits) {
        if (numTraits > 1) {
            throw new RuntimeException("not yet implemented");
        }

        double[] belowTraits = new double[nObservedTips * dimTrait];
        int nTaxa = getParameter().getParameterCount();

        DenseMatrix64F P = DenseMatrix64F.wrap(dimTrait, dimTrait, samplingPrecisionParameter.getParameterValues());
        DenseMatrix64F Q = new DenseMatrix64F(dimTrait, dimTrait);
        DenseMatrix64F V = new DenseMatrix64F(dimTrait, dimTrait);

        DenseMatrix64F P0 = new DenseMatrix64F(dimTrait, dimTrait);

        int[] wrappedIndices = new int[dimTrait];
        for (int i = 0; i < dimTrait; i++) {
            wrappedIndices[i] = i;
        }

        WrappedVector n = new WrappedVector.Raw(new double[dimTrait]);

        int precisionOffset = precisionType.getPrecisionOffset(dimTrait);
        int meanOffset = precisionType.getMeanOffset(dimTrait);
        int repOffset = precisionType.getPartialsDimension(dimTrait);
        int dimPrecision = precisionType.getPrecisionLength(dimTrait);

        int aboveOffset = 0;
        int belowOffset = 0;
        for (int i = 0; i < nTaxa; i++) {
            double[] partial = childModel.getTipPartial(i, false);
            WrappedVector.Indexed x = new WrappedVector.Indexed(aboveTraits, aboveOffset, wrappedIndices);

            for (int j : relevantRepeats[i]) {
                System.arraycopy(partial, j * repOffset + precisionOffset, P0.data, 0, dimPrecision);
                WrappedVector.Indexed m0 = new WrappedVector.Indexed(partial, j * repOffset + meanOffset, wrappedIndices);


                boolean completelyObserved = true;
                for (int k = 0; k < dimTrait; k++) {
                    if (P0.get(k, k) < Double.POSITIVE_INFINITY) {
                        completelyObserved = false;
                        break;
                    }
                }

                if (completelyObserved) {
                    for (int k = 0; k < dimTrait; k++) {
                        belowTraits[belowOffset + k] = m0.get(k);
                    }
                } else {
                    CommonOps.add(P0, P, Q);
                    MissingOps.safeInvert2(Q, V, false);

                    MissingOps.safeWeightedAverage(m0, P0, x, P, n, V, dimTrait);

                    double[] sample = MissingOps.nextPossiblyDegenerateNormal(n, V);

                    System.arraycopy(sample, 0, belowTraits, belowOffset, dimTrait);
                }

                belowOffset += dimTrait;

            }

            aboveOffset += dimTrait;


        }

        return belowTraits;
    }

    @Override
    public double[] transformTreeTraits(double[] treeTraits) {
        double[] repeatedTraits = new double[dimTrait * nObservedTips];
        int originalOffset = 0;
        int expandedOffset = 0;
        for (ArrayList<Integer> repeats : relevantRepeats) {
            for (int i : repeats) {
                System.arraycopy(treeTraits, originalOffset, repeatedTraits, expandedOffset, dimTrait);
                expandedOffset += dimTrait;
            }
            originalOffset += dimTrait;
        }

        return repeatedTraits;
    }

    private static final boolean DEBUG = false;

    // TODO Move remainder into separate class file
    public static final String REPEATED_MEASURES_MODEL = "repeatedMeasuresModel";
    private static final String PRECISION = "samplingPrecision";
    private static final String SCALE_BY_TIP_HEIGHT = "scaleByTipHeight";

    public static AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {
        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            final ContinuousTraitPartialsProvider subModel;


            XMLObject cxo = xo.getChild(PRECISION);
            MatrixParameterInterface samplingPrecision = (MatrixParameterInterface)
                    cxo.getChild(MatrixParameterInterface.class);

            CholeskyDecomposition chol;
            try {
                chol = new CholeskyDecomposition(samplingPrecision.getParameterAsMatrix());
            } catch (IllegalDimension illegalDimension) {
                throw new XMLParseException(PRECISION + " must be a square matrix.");
            }

            if (!chol.isSPD()) {
                throw new XMLParseException(PRECISION + " must be a positive definite matrix.");
            }


            boolean scaleByTipHeight = xo.getAttribute(SCALE_BY_TIP_HEIGHT, false);

            int dimTrait = samplingPrecision.getColumnDimension();
            final PrecisionType precisionType;
            if (xo.getAttribute(ContinuousTraitDataModelParser.FORCE_FULL_PRECISION, false) ||
                    dimTrait > 1) {
                precisionType = PrecisionType.FULL;
            } else {
                precisionType = PrecisionType.SCALAR;
            }

            if (xo.hasChildNamed(TreeTraitParserUtilities.TRAIT_PARAMETER)) {
                subModel = ContinuousTraitDataModelParser.parseContinuousTraitDataModel(xo, precisionType);
            } else {
                subModel = (ContinuousTraitPartialsProvider) xo.getChild(ContinuousTraitPartialsProvider.class);
                if (subModel.getPrecisionType() != precisionType) {
                    throw new XMLParseException("Precision type of " + REPEATED_MEASURES_MODEL + " is " +
                            precisionType.getClass() + ", but the precision type of the child model " +
                            subModel.getModelName() + " is " + subModel.getPrecisionType().getClass());
                }
            }
            String modelName = subModel.getModelName();

            int numTraits = xo.getAttribute(NUM_TRAITS, subModel.getTraitCount());

            if (subModel.getTraitDimension() != dimTrait) {
                throw new XMLParseException("sub-model has trait dimension " + subModel.getTraitDimension() +
                        ", but sampling precision has dimension " + dimTrait);
            }

            if (!scaleByTipHeight) {
                return new RepeatedMeasuresTraitDataModel(
                        modelName,
                        subModel,
                        subModel.getParameter(),
                        subModel.getDataMissingIndicators(),
//                    missingIndicators,
                        true,
                        dimTrait,
                        numTraits,
//                    diffusionModel.getPrecisionParameter().getRowDimension(),
                        samplingPrecision,
                        precisionType
                );
            } else {
                return new TreeScaledRepeatedMeasuresTraitDataModel(
                        modelName,
                        subModel,
                        subModel.getParameter(),
                        subModel.getDataMissingIndicators(),
                        true,
                        dimTrait,
                        subModel.getTraitCount(),
                        samplingPrecision,
                        precisionType
                );
            }
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
            return RepeatedMeasuresTraitDataModel.class;
        }

        @Override
        public String getParserName() {
            return REPEATED_MEASURES_MODEL;
        }
    };

    private final static XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(PRECISION, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class),
            }),
            // Tree trait parser
//            new ElementRule(MutableTreeModel.class),
//            AttributeRule.newStringRule(TreeTraitParserUtilities.TRAIT_NAME),
            new XORRule(
                    new ElementRule(ContinuousTraitPartialsProvider.class),
                    new AndRule(ContinuousTraitDataModelParser.rules)
            ),
            new ElementRule(TreeTraitParserUtilities.MISSING, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }, true),
            AttributeRule.newBooleanRule(SCALE_BY_TIP_HEIGHT, true),
//            new ElementRule(MultivariateDiffusionModel.class),
    };


}
