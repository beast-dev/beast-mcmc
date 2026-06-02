/*
 * RepeatedMeasuresTraitDataModelParser.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.evomodelxml.treedatalikelihood.continuous;

import dr.app.beauti.components.BeautiModelIDProvider;
import dr.app.beauti.components.BeautiParameterIDProvider;
import dr.evomodel.treedatalikelihood.continuous.ContinuousTraitPartialsProvider;
import dr.evomodel.treedatalikelihood.continuous.RepeatedMeasuresTraitDataModel;
import dr.evomodel.treedatalikelihood.continuous.TreeScaledRepeatedMeasuresTraitDataModel;
import dr.evomodel.treedatalikelihood.continuous.cdi.PrecisionType;
import dr.evomodelxml.continuous.ContinuousTraitDataModelParser;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;
import dr.math.matrixAlgebra.CholeskyDecomposition;
import dr.math.matrixAlgebra.IllegalDimension;
import dr.xml.*;

import static dr.evomodelxml.continuous.ContinuousTraitDataModelParser.NUM_TRAITS;
import static dr.evomodelxml.treedatalikelihood.ContinuousDataLikelihoodParser.FORCE_FULL_PRECISION;


public class RepeatedMeasuresTraitDataModelParser extends AbstractXMLObjectParser implements BeautiModelIDProvider {
    public static final String REPEATED_MEASURES_MODEL = "repeatedMeasuresModel";
    private static final String PRECISION = "samplingPrecision";
    private static final String SCALE_BY_TIP_HEIGHT = "scaleByTipHeight";


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

        // Jitter
        TreeTraitParserUtilities utilities = new TreeTraitParserUtilities(); // TODO: ideally this wouldn't be here
        if (xo.hasChildNamed(TreeTraitParserUtilities.JITTER)) {
            utilities.jitter(xo, samplingPrecision.getColumnDimension(), subModel.getDataMissingIndicators());
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
            TreeTraitParserUtilities.jitterRules(true),
            AttributeRule.newBooleanRule(FORCE_FULL_PRECISION, true),
    };

    //********************************************************************
    // BeautiModelIDProvider interface
    //********************************************************************

    public String getParserTag() {
        return REPEATED_MEASURES_MODEL;
    }

    public String getId(String modelName) {
        return modelName + ".residualModel";
    }

//    public String getId(ContinuousModelExtensionType extensionType, String modelName) {
//        String extendedName;
//        switch (extensionType) {
//            case RESIDUAL:
//                extendedName = "residualModel";
//                break;
//            case LATENT_FACTORS:
//                extendedName = "factorModel";
//                break;
//            case NONE:
//                throw new IllegalArgumentException("Should not be called");
//            default:
//                throw new IllegalArgumentException("Unknown extension type");
//
//        }

//        return modelName + "." + extendedName;
//    }

    private static final String PRECISION_ID = "samplingPrecision";

    public static final String EXTENSION_PRECISION = "extensionPrecision";
    public static final String EXTENSION_VARIANCE = "extensionVarCovar";

    private static final String[] ALLOWABLE_PARAMETERS = new String[]{EXTENSION_PRECISION, EXTENSION_VARIANCE};

//    private final BeautiParameterIDProvider extensionPrecisionIDProvider = new BeautiParameterIDProvider(EXTENSION_PRECISION);
//    private final BeautiParameterIDProvider extensionPrecisionIDProvider = new BeautiParameterIDProvider(EXTENSION_PRECISION);


    public BeautiParameterIDProvider getBeautiParameterIDProvider(String parameterKey) {
        for (int i = 0; i < ALLOWABLE_PARAMETERS.length; i++) {
            if (parameterKey.equalsIgnoreCase(ALLOWABLE_PARAMETERS[i])) {
                return new BeautiParameterIDProvider(parameterKey);
            }
        }
        throw new IllegalArgumentException("Unrecognized parameter key '" + parameterKey + "'");
    }
}
