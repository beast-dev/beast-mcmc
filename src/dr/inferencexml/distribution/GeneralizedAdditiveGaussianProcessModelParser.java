/*
 * GeneralizedAdditiveGaussianProcessModelParser.java
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

package dr.inferencexml.distribution;

import dr.inference.distribution.*;
import dr.inference.model.DesignMatrix;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

//import static dr.inferencexml.distribution.GeneralizedLinearModelParser.*;
//import static dr.inferencexml.glm.ExperimentalGeneralizedLinearModelParser.DEPENDENT_VARIABLES;

/**
 * @author Filippo Monti
 * @author Marc Suchard
 */
public class GeneralizedAdditiveGaussianProcessModelParser extends AbstractXMLObjectParser {

    public static final String GAM_GP_LIKELIHOOD = "gamGpModel";
    private static final String REALIZED_FIELD = "realizedField";
    private static final String GAUSSIAN_NOISE = "gaussianNoise";
    private static final String KERNEL = "kernel";
    public static final String DEPENDENT_VARIABLES = "dependentVariables";
    public static final String INDEPENDENT_VARIABLES = "independentVariables";
    public static final String INDICATOR = "indicator";

    public String getParserName() {
        return GAM_GP_LIKELIHOOD;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

//        XMLObject cxo = xo.getChild(DEPENDENT_VARIABLES);
//        Parameter dependentParam = null;
//        if (cxo != null) {
//            dependentParam = (Parameter) cxo.getChild(Parameter.class);
//        }

//        Parameter realizedField = (Parameter) xo.getElementFirstChild(REALIZED_FIELD);

        Parameter realizedField = (Parameter) xo.getElementFirstChild(REALIZED_FIELD);
        Parameter gaussianNoise = (Parameter) xo.getElementFirstChild(GAUSSIAN_NOISE);
        Parameter kernelParameter = (Parameter) xo.getElementFirstChild(KERNEL);

        List<DesignMatrix> matrices = new ArrayList<>();
        for (XMLObject cxo : xo.getAllChildren(INDEPENDENT_VARIABLES)) {
            DesignMatrix designMatrix = (DesignMatrix) cxo.getChild(DesignMatrix.class);
            matrices.add(designMatrix);
        }

//        String family = xo.getStringAttribute(FAMILY);
        LogGaussianProcessModel gp = new LogGaussianProcessModel(realizedField, gaussianNoise, kernelParameter, matrices);

        addIndependentParameters(xo, gp, null);

        return gp;
    }

//    public void addRandomEffects(XMLObject xo, GeneralizedLinearModel glm,
//                                 Parameter dependentParam) throws XMLParseException {
//        int totalCount = xo.getChildCount();
//
//        for (int i = 0; i < totalCount; i++) {
//            if (xo.getChildName(i).compareTo(RANDOM_EFFECTS) == 0) {
//                XMLObject cxo = (XMLObject) xo.getChild(i);
//                Parameter randomEffect = (Parameter) cxo.getChild(Parameter.class);
//                checkRandomEffectsDimensions(randomEffect, dependentParam);
//                glm.addRandomEffectsParameter(randomEffect);
//            }
//        }
//    }

    // TODO remove code duplication
    public void addIndependentParameters(XMLObject xo, LogGaussianProcessModel lgpm,
                                         Parameter dependentParam) throws XMLParseException {
        int totalCount = xo.getChildCount();

        for (XMLObject cxo : xo.getAllChildren(INDEPENDENT_VARIABLES)) {
//            Parameter field = (Parameter) cxo.getChild(Parameter.class);
            DesignMatrix designMatrix = (DesignMatrix) cxo.getChild(DesignMatrix.class);

//            if (field.getDimension() != designMatrix.getRowDimension()) {
//                throw new XMLParseException(field.getDimension() + " " + designMatrix.getRowDimension());
//            }

            lgpm.addIndependentParameter(designMatrix, designMatrix, null);
//            glm.addIndependentParameter(field, designMatrix, null);

        }

//        for (int i = 0; i < totalCount; i++) {
//            if (xo.getChildName(i).compareTo(INDEPENDENT_VARIABLES) == 0) {
//                XMLObject cxo = (XMLObject) xo.getChild(i);
//                Parameter independentParam = (Parameter) cxo.getChild(Parameter.class);
//                DesignMatrix designMatrix = (DesignMatrix) cxo.getChild(DesignMatrix.class);
//                checkDimensions(independentParam, dependentParam, designMatrix);
//                cxo = cxo.getChild(INDICATOR);
//                Parameter indicator = null;
//                if (cxo != null) {
//                    indicator = (Parameter) cxo.getChild(Parameter.class);
//                    if (indicator.getDimension() <= 1) {
//                        // if a dimension hasn't been set, then set it automatically
//                        indicator.setDimension(independentParam.getDimension());
//                    }
//                    if (indicator.getDimension() != independentParam.getDimension())
//                        throw new XMLParseException("dim(" + independentParam.getId() + ") != dim(" + indicator.getId() + ")");
//                }
//
////                if (checkFullRankOfMatrix) {
////                    checkFullRank(designMatrix);
////                }
//
//                glm.addIndependentParameter(independentParam, designMatrix, indicator);
//            }
//        }
    }

    // TODO remove code duplication
    private void checkDimensions(Parameter independentParam, Parameter dependentParam, DesignMatrix designMatrix)
            throws XMLParseException {
        if (dependentParam != null) {
            if (dependentParam.getDimension() <= 1) {
                // if a dimension hasn't been set, then set it automatically
                dependentParam.setDimension(designMatrix.getRowDimension());
            }
            if ((dependentParam.getDimension() != designMatrix.getRowDimension()) ||
                    (independentParam.getDimension() != designMatrix.getColumnDimension()))
                throw new XMLParseException(
                        "dim(" + dependentParam.getId() + ") != dim(" + designMatrix.getId() + " %*% " + independentParam.getId() + ")"
                );
        } else {
            if (independentParam.getDimension() <= 1) {
                // if a dimension hasn't been set, then set it automatically
                independentParam.setDimension(designMatrix.getColumnDimension());
            }
            if (independentParam.getDimension() != designMatrix.getColumnDimension()) {
                throw new XMLParseException(
                        "dim(" + independentParam.getId() + ") is incompatible with dim (" + designMatrix.getId() + ")"
                );
            }
//            System.err.println(independentParam.getId()+" and "+designMatrix.getId());
        }
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
//            AttributeRule.newStringRule(FAMILY),
//            AttributeRule.newBooleanRule(CHECK_IDENTIFIABILITY, true),
//            AttributeRule.newBooleanRule(CHECK_FULL_RANK, true),
            new ElementRule(REALIZED_FIELD,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }),
            new ElementRule(GAUSSIAN_NOISE,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }),
            new ElementRule(KERNEL,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }),
            new ElementRule(DEPENDENT_VARIABLES,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
//            new ElementRule(REALIZED_FIELD,
//                    new XMLSyntaxRule[] { new ElementRule(Parameter.class) } ),
            new ElementRule(INDEPENDENT_VARIABLES,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class, true),
                            new ElementRule(DesignMatrix.class),
                            new ElementRule(INDICATOR,
                                    new XMLSyntaxRule[]{
                                            new ElementRule(Parameter.class)
                                    }, true),
                    }, 0, 10),
//            new ElementRule(RANDOM_EFFECTS,
//                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, 0, 3),
    };

    public String getParserDescription() {
        return "Calculates the generalized linear model likelihood of the dependent parameters given one or more blocks of independent parameters and their design matrix.";
    }

    public Class getReturnType() {
        return LogGaussianProcessModel.class;
    }
}
