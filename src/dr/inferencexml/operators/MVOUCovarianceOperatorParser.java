/*
 * MVOUCovarianceOperatorParser.java
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

package dr.inferencexml.operators;

import dr.inference.model.MatrixParameter;
import dr.inference.operators.AdaptableMCMCOperator;
import dr.inference.operators.AdaptationMode;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.MVOUCovarianceOperator;
import dr.xml.*;

/**
 *
 */
public class MVOUCovarianceOperatorParser extends AbstractXMLObjectParser {

    public static final String MVOU_OPERATOR = "mvouOperator";
    public static final String MIXING_FACTOR = "mixingFactor";
    public static final String VARIANCE_MATRIX = "varMatrix";
    public static final String PRIOR_DF = "priorDf";

    public String getParserName() {
        return MVOU_OPERATOR;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        AdaptationMode mode = AdaptationMode.parseMode(xo);
        double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        double mixingFactor = xo.getDoubleAttribute(MIXING_FACTOR);
        int priorDf = xo.getIntegerAttribute(PRIOR_DF);

        if (mixingFactor <= 0.0 || mixingFactor >= 1.0) {
            throw new XMLParseException("mixingFactor must be greater than 0.0 and less thatn 1.0");
        }

//            Parameter parameter = (Parameter) xo.getChild(Parameter.class);

//            XMLObject cxo = (XMLObject) xo.getChild(VARIANCE_MATRIX);
        MatrixParameter varMatrix = (MatrixParameter) xo.getChild(MatrixParameter.class);

        // Make sure varMatrix is square and dim(varMatrix) = dim(parameter)

        if (varMatrix.getColumnDimension() != varMatrix.getRowDimension())
            throw new XMLParseException("The variance matrix is not square");

//            if (varMatrix.getColumnDimension() != parameter.getDimension())
//                throw new XMLParseException("The parameter and variance matrix have differing dimensions");

        return new MVOUCovarianceOperator(mixingFactor, varMatrix, priorDf, weight, mode);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element returns junk.";
    }

    public Class getReturnType() {
        return MVOUCovarianceOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            AttributeRule.newDoubleRule(MIXING_FACTOR),
            AttributeRule.newIntegerRule(PRIOR_DF),
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            AttributeRule.newBooleanRule(AdaptableMCMCOperator.AUTO_OPTIMIZE, true),
//                new ElementRule(Parameter.class),
//                new ElementRule(VARIANCE_MATRIX,
//                        new XMLSyntaxRule[]{new ElementRule(MatrixParameter.class)}),

            new ElementRule(MatrixParameter.class)

    };

}
