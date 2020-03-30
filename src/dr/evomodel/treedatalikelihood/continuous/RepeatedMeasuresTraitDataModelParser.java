/*
 * RepeatedMeasuresTraitDataModelParser.java
 *
 * Copyright (c) 2002-2020 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.app.beauti.components.BeautiModelIDProvider;
import dr.app.beauti.components.BeautiParameterIDProvider;
import dr.app.beauti.components.continuous.ContinuousModelExtensionType;
import dr.evolution.tree.MutableTreeModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.inference.model.CompoundParameter;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;
import dr.math.matrixAlgebra.CholeskyDecomposition;
import dr.math.matrixAlgebra.IllegalDimension;
import dr.xml.*;

import java.util.List;

public class RepeatedMeasuresTraitDataModelParser extends AbstractXMLObjectParser implements BeautiModelIDProvider {
    public static final String REPEATED_MEASURES_MODEL = "repeatedMeasuresModel";
    private static final String PRECISION = "samplingPrecision";

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

    //********************************************************************
    // BeautiModelIDProvider interface
    //********************************************************************

    public String getParserTag() {
        return REPEATED_MEASURES_MODEL;
    }

    public String getId(String modelName) {
        throw new IllegalArgumentException("Should not be called");
    }

    public String getId(ContinuousModelExtensionType extensionType, String modelName) {
        String extendedName;
        switch (extensionType) {
            case RESIDUAL:
                extendedName = "residualModel";
                break;
            case LATENT_FACTORS:
                extendedName = "factorModel";
                break;
            case NONE:
                throw new IllegalArgumentException("Should not be called");
            default:
                throw new IllegalArgumentException("Unknown extension type");

        }

        return modelName + "." + extendedName;
    }

    private static final String PRECISION_ID = "samplingPrecision";

    private final BeautiParameterIDProvider extensionPrecisionIDProvider = new BeautiParameterIDProvider("extensionPrecision");

    public BeautiParameterIDProvider getBeautiParameterIDProvider(String parameterKey) {
        assert parameterKey.equals("extensionPrecision")
                : "Only the 'extensionPrecision' parameter is implemented for the 'repeatedMeasuresModel'.";

        return extensionPrecisionIDProvider;
    }
}
