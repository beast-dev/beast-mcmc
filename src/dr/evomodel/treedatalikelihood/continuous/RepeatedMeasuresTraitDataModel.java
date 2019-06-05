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

import dr.evolution.tree.MutableTreeModel;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.continuous.cdi.PrecisionType;
import dr.evomodel.treedatalikelihood.preorder.ContinuousExtensionDelegate;
import dr.evomodel.treedatalikelihood.preorder.ModelExtensionProvider;
import dr.evomodel.treedatalikelihood.preorder.ProcessSimulationDelegate;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.inference.model.CompoundParameter;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.math.matrixAlgebra.CholeskyDecomposition;
import dr.math.matrixAlgebra.IllegalDimension;
import dr.math.matrixAlgebra.Matrix;
import dr.math.matrixAlgebra.missingData.MissingOps;
import dr.xml.*;
import org.ejml.data.DenseMatrix64F;

import java.util.List;

/**
 * @author Marc A. Suchard
 * @author Gabriel Hassler
 */
public class RepeatedMeasuresTraitDataModel extends
        ContinuousTraitDataModel implements ContinuousTraitPartialsProvider, ModelExtensionProvider.NormalExtensionProvider {

    private final String traitName;
    private final MatrixParameterInterface samplingPrecision;
    private boolean diagonalOnly = false;
    private Matrix samplingVariance;
    private boolean varianceKnown = false;

    public RepeatedMeasuresTraitDataModel(String name,
                                          CompoundParameter parameter,
                                          List<Integer> missingIndices,
//                                          boolean[] missindIndicators,
                                          boolean useMissingIndices,
                                          final int dimTrait,
                                          MatrixParameterInterface samplingPrecision) {
        super(name, parameter, missingIndices, useMissingIndices, dimTrait, PrecisionType.FULL);
        this.traitName = name;
        this.samplingPrecision = samplingPrecision;
        addVariable(samplingPrecision);

//        this.samplingVariance = new Matrix(samplingPrecision.getParameterAsMatrix()).inverse();
        this.samplingVariance = null;


        samplingPrecision.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0,
                samplingPrecision.getDimension()));

    }

    @Override
    public double[] getTipPartial(int taxonIndex, boolean fullyObserved) {

        assert (numTraits == 1);
        assert (samplingPrecision.getRowDimension() == dimTrait && samplingPrecision.getColumnDimension() == dimTrait);

        recomputeVariance();

        if (fullyObserved) {
            return new double[dimTrait + 1];
        }

        double[] partial = super.getTipPartial(taxonIndex, fullyObserved);
        DenseMatrix64F V = MissingOps.wrap(partial, dimTrait + dimTrait * dimTrait, dimTrait, dimTrait);

        //TODO: remove diagonalOnly part
        if (diagonalOnly) {
            for (int index = 0; index < dimTrait; index++) {
                V.set(index, index, V.get(index, index) + 1 / samplingPrecision.getParameterValue(index));
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

        return partial;
    }

    public MatrixParameterInterface getPrecisionMatrix() {
        return samplingPrecision;
    }


    public Parameter getSamplingPrecision() {
        return samplingPrecision;
    }

    private void recomputeVariance() {
        if (!varianceKnown) {
            samplingVariance = new Matrix(samplingPrecision.getParameterAsMatrix()).inverse();
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

        if (variable == samplingPrecision) {

            varianceKnown = false;
            fireModelChanged();
        }
    }

    @Override
    public ContinuousExtensionDelegate getExtensionDelegate(ContinuousDataLikelihoodDelegate delegate, TreeTrait treeTrait, Tree tree) {
        return new ContinuousExtensionDelegate.MultivariateNormalExtensionDelegate(delegate, treeTrait, this, tree);
    }

    @Override
    public DenseMatrix64F getExtensionVariance() {
        double[] buffer = samplingVariance.toArrayComponents();
        return DenseMatrix64F.wrap(dimTrait, dimTrait, buffer);
    }

    @Override
    public boolean[] getMissingVector() {
        return super.getMissingIndicator();
    }

    // TODO Move remainder into separate class file
    private static final String REPEATED_MEASURES_MODEL = "repeatedMeasuresModel";
    private static final String PRECISION = "samplingPrecision";

    public static AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {
        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            MutableTreeModel treeModel = (MutableTreeModel) xo.getChild(TreeModel.class);
            TreeTraitParserUtilities utilities = new TreeTraitParserUtilities();

            TreeTraitParserUtilities.TraitsAndMissingIndices returnValue =
                    utilities.parseTraitsFromTaxonAttributes(xo, TreeTraitParserUtilities.DEFAULT_TRAIT_NAME,
                            treeModel, true);
            CompoundParameter traitParameter = returnValue.traitParameter;
            List<Integer> missingIndices = returnValue.missingIndices;

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


            String traitName = returnValue.traitName;
            //TODO diffusionModel was only used for the dimension.
            // But this should be the same as the samplingPrecision dimension ?
//            MultivariateDiffusionModel diffusionModel = (MultivariateDiffusionModel)
//                    xo.getChild(MultivariateDiffusionModel.class);

            //TODO: This was never used.
//            final boolean[] missingIndicators = new boolean[returnValue.traitParameter.getDimension()];
//            for (int i : missingIndices) {
//                missingIndicators[i] = true;
//            }

            return new RepeatedMeasuresTraitDataModel(
                    traitName,
                    traitParameter,
                    missingIndices,
//                    missingIndicators,
                    true,
                    samplingPrecision.getColumnDimension(),
//                    diffusionModel.getPrecisionParameter().getRowDimension(),
                    samplingPrecision
            );
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
            new ElementRule(MutableTreeModel.class),
            AttributeRule.newStringRule(TreeTraitParserUtilities.TRAIT_NAME),
            new ElementRule(TreeTraitParserUtilities.TRAIT_PARAMETER, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }),
            new ElementRule(TreeTraitParserUtilities.MISSING, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }, true),
//            new ElementRule(MultivariateDiffusionModel.class),
    };


}
