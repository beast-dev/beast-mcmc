/*
 * CorrelationDiagonalMatrixParser.java
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

package dr.inferencexml.model;

import dr.inference.model.*;
import dr.xml.*;

/**
 *
 */
public class DiagonalContrainedMatrixViewParser extends AbstractXMLObjectParser {

    public final static String VIEW_NAME = "diagonalContrainedMatrixView";
    public static final String CONSTRAINT_VALUE = "constraintValue";

    public String getParserName() {
        return VIEW_NAME;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter mask = (Parameter) xo.getElementFirstChild(MaskedParameterParser.MASKING);
        double constraintValue = xo.getAttribute(CONSTRAINT_VALUE, 1.0);
        MatrixParameterInterface matrix = (MatrixParameterInterface)
                xo.getChild(MatrixParameterInterface.class);

        if (matrix.getRowDimension() != matrix.getColumnDimension()
                && matrix.getRowDimension() != mask.getDimension()) {
            throw new XMLParseException("Incompatiable matrix and mask dimensions");
        }

        return new DiagonalConstrainedMatrixView(xo.getId(), matrix, mask, constraintValue);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A diagonal matrix parameter constructed from its diagonals.";
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            AttributeRule.newDoubleRule(CONSTRAINT_VALUE, true),
            new ElementRule(MatrixParameterInterface.class),
            new ElementRule(MaskedParameterParser.MASKING,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }),
    };

    public Class getReturnType() {
        return MatrixParameter.class;
    }
    
}
