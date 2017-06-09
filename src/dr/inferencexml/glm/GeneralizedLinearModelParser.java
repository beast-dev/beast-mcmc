/*
 * GeneralizedLinearModelParser.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.inferencexml.glm;

import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.SingularValueDecomposition;
import dr.inference.distribution.DensityModel;
import dr.inference.distribution.ParametricDistributionModel;
import dr.inference.glm.*;
import dr.inference.model.DesignMatrix;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.util.Transform;
import dr.xml.*;

/**
 *
 */
public class GeneralizedLinearModelParser extends AbstractXMLObjectParser {

    public static final String GLM_LIKELIHOOD = "glmModelNew";

    public static final String DEPENDENT_VARIABLES = "dependentVariables";
    public static final String INDEPENDENT_VARIABLES = "independentVariables";

    public static final String MODEL = "model";
    public static final String LINK_FUNCTION = "linkFunction";

    public static final String INDICATOR = "indicator";
    public static final String RANDOM_EFFECTS = "randomEffects";
    public static final String CHECK_IDENTIFIABILITY = "checkIdentifiability";
    public static final String CHECK_FULL_RANK = "checkFullRank";

    public String getParserName() {
        return GLM_LIKELIHOOD;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {


        XMLObject cxo = xo.getChild(MODEL);
        String linkFunctionName = cxo.getAttribute(LINK_FUNCTION, "identity");
        DensityModel model = (DensityModel) cxo.getChild(DensityModel.class);
        GeneralizedLinearModel.LinkFunction linkFunction = GeneralizedLinearModel.LinkFunction.valueOf(linkFunctionName.toUpperCase());

        cxo = xo.getChild(DEPENDENT_VARIABLES);
        Parameter dependentParam = (Parameter) cxo.getChild(Parameter.class);

        GeneralizedLinearModel glm = new GeneralizedLinearModel(dependentParam, model, linkFunction);

        // parse this as parameters in a DistributionModel?
//        if (glm.requiresScale()) {
//            cxo = xo.getChild(SCALE_VARIABLES);
//            Parameter scaleParameter = null;
////                DesignMatrix designMatrix = null;
//            Parameter scaleDesign = null;
//            if (cxo != null) {
//                scaleParameter = (Parameter) cxo.getChild(Parameter.class);
//                XMLObject gxo = cxo.getChild(INDICATOR);
//                if (gxo != null)
//                    scaleDesign = (Parameter) gxo.getChild(Parameter.class);
////                    designMatrix = (DesignMatrix) cxo.getChild(DesignMatrix.class);
//            }
//            if (scaleParameter == null)
//                throw new XMLParseException("Family '" + family + "' requires scale parameters");
//            if (scaleDesign == null)
//                scaleDesign = new Parameter.Default(dependentParam.getDimension(), 0.0);
//            else {
//                if (scaleDesign.getDimension() != dependentParam.getDimension())
//                    throw new XMLParseException("Scale ("+dependentParam.getDimension()+") and scaleDesign parameters ("+scaleDesign.getDimension()+") must be the same dimension");
//                for (int i = 0; i < scaleDesign.getDimension(); i++) {
//                    double value = scaleDesign.getParameterValue(i);
//                    if (value < 1 || value > scaleParameter.getDimension())
//                        throw new XMLParseException("Invalid scaleDesign value");
//                    scaleDesign.setParameterValue(i, value - 1);
//                }
//            }
//
//            glm.addScaleParameter(scaleParameter, scaleDesign);
//        }


//        System.err.println("START 0");
        addIndependentParameters(xo, glm, dependentParam);
//        System.err.println("START 1");
        addRandomEffects(xo, glm, dependentParam);
//        System.err.println("START 2");

        boolean checkIdentifiability = xo.getAttribute(CHECK_IDENTIFIABILITY, true);
        if (checkIdentifiability) {
            if (!glm.getAllIndependentVariablesIdentifiable()) {
                throw new XMLParseException("All design matrix predictors are not identifiable in "+  xo.getId());
            }
        }
//        System.err.println("PASSED B");
        checkFullRankOfMatrix = xo.getAttribute(CHECK_FULL_RANK,true);
//        System.err.println("PASSED C");
        return glm;
    }

    public void addRandomEffects(XMLObject xo, GeneralizedLinearModel glm,
                                 Parameter dependentParam) throws XMLParseException {
        int totalCount = xo.getChildCount();

        for (int i = 0; i < totalCount; i++) {
            if (xo.getChildName(i).compareTo(RANDOM_EFFECTS) == 0) {
                XMLObject cxo = (XMLObject) xo.getChild(i);
                Parameter randomEffect = (Parameter) cxo.getChild(Parameter.class);
                checkRandomEffectsDimensions(randomEffect, dependentParam);
                glm.addRandomEffectsParameter(randomEffect);
            }
        }
    }

    public void addIndependentParameters(XMLObject xo, GeneralizedLinearModel glm,
                                         Parameter dependentParam) throws XMLParseException {
        int totalCount = xo.getChildCount();
//        System.err.println("number of independent parameters = "+totalCount);

        for (int i = 0; i < totalCount; i++) {
            if (xo.getChildName(i).compareTo(INDEPENDENT_VARIABLES) == 0) {
                XMLObject cxo = (XMLObject) xo.getChild(i);
                Parameter independentParam = (Parameter) cxo.getChild(Parameter.class);
                DesignMatrix designMatrix = (DesignMatrix) cxo.getChild(DesignMatrix.class);
                checkDimensions(independentParam, dependentParam, designMatrix);
                cxo = cxo.getChild(INDICATOR);
                Parameter indicator = null;
                if (cxo != null) {
                    indicator = (Parameter) cxo.getChild(Parameter.class);
                    if (indicator.getDimension() != independentParam.getDimension())
                        throw new XMLParseException("dim(" + independentParam.getId() + ") != dim(" + indicator.getId() + ")");
                }
//                System.err.println("A");
                if (checkFullRankOfMatrix) {
                    checkFullRank(designMatrix);
                }
//                System.err.println("B");

//                System.err.println(new Matrix(designMatrix.getParameterAsMatrix()));
//                System.exit(-1);

                glm.addIndependentParameter(independentParam, designMatrix, indicator);
//                System.err.println("C");
            }
        }
    }

    private boolean checkFullRankOfMatrix;

    private void checkFullRank(DesignMatrix designMatrix) throws XMLParseException {
        int fullRank = designMatrix.getColumnDimension();
//        System.err.println("designMatrix getColumnDimension = "+fullRank);
        SingularValueDecomposition svd = new SingularValueDecomposition(
                new DenseDoubleMatrix2D(designMatrix.getParameterAsMatrix()));
        int realRank = svd.rank();
        if (realRank != fullRank) {
            throw new XMLParseException(
                    "rank(" + designMatrix.getId() + ") = " + realRank +
                            ".\nMatrix is not of full rank as colDim(" + designMatrix.getId() + ") = " + fullRank
            );
        }
    }

    private void checkRandomEffectsDimensions(Parameter randomEffect, Parameter dependentParam)
            throws XMLParseException {
        if (dependentParam != null) {
            if (randomEffect.getDimension() != dependentParam.getDimension()) {
                throw new XMLParseException(
                        "dim(" + dependentParam.getId() + ") != dim(" + randomEffect.getId() + ")"
                );
            }
        }
    }

    private void checkDimensions(Parameter independentParam, Parameter dependentParam, DesignMatrix designMatrix)
            throws XMLParseException {
        if (dependentParam != null) {
            if ((dependentParam.getDimension() != designMatrix.getRowDimension()) ||
                    (independentParam.getDimension() != designMatrix.getColumnDimension()))
                throw new XMLParseException(
                        "dim(" + dependentParam.getId() + ") != dim(" + designMatrix.getId() + " %*% " + independentParam.getId() + ")"
                );
        } else {
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
            AttributeRule.newBooleanRule(CHECK_IDENTIFIABILITY, true),
            AttributeRule.newBooleanRule(CHECK_FULL_RANK, true),
            new ElementRule(MODEL,
                    new XMLSyntaxRule[]{
                            AttributeRule.newStringRule(LINK_FUNCTION),
                            new ElementRule(ParametricDistributionModel.class),
                    }),
            new ElementRule(DEPENDENT_VARIABLES,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, false),
            new ElementRule(INDEPENDENT_VARIABLES,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class, true),
                            new ElementRule(DesignMatrix.class),
                            new ElementRule(INDICATOR,
                                    new XMLSyntaxRule[]{
                                            new ElementRule(Parameter.class)
                                    }, true),
                    }, 1, 3),
            new ElementRule(RANDOM_EFFECTS,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, 0, 3),
//				new ElementRule(BASIS_MATRIX,
//						new XMLSyntaxRule[]{new ElementRule(DesignMatrix.class)})
    };

    public String getParserDescription() {
        return "Calculates the generalized linear model likelihood of the dependent parameters given one or more blocks of independent parameters and their design matrix.";
    }

    public Class getReturnType() {
        return GeneralizedLinearModel.class;
    }
}
