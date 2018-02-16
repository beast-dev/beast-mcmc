/*
 * CompoundSymmetricMatrixParser.java
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

package dr.inferencexml.model;

import dr.inference.model.CorrelationSymmetricMatrix;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 *
 */
public class CorrelationSymmetricMatrixParser extends AbstractXMLObjectParser {

    private final static String MATRIX_PARAMETER = "compoundSymmetricMatrix";
    private static final String DIAGONAL = "diagonal";
    private static final String OFF_DIAGONAL = "offDiagonal";

    public String getParserName() {
        return MATRIX_PARAMETER;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter diagonals = (Parameter) xo.getElementFirstChild(DIAGONAL);
        Parameter offDiagonals = (Parameter) xo.getElementFirstChild(OFF_DIAGONAL);

        if (offDiagonals.getDimension() !=
                (diagonals.getDimension() * (diagonals.getDimension() - 1) / 2)) {
            throw new XMLParseException("Invalid parameter dimensions in `" + xo.getId() + "'");
        }

        return new CorrelationSymmetricMatrix(diagonals, offDiagonals);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A (partial) correlation-parameterized matrix parameter constructed from its diagonals and off-diagonals.";
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(DIAGONAL,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(OFF_DIAGONAL,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
    };

    public Class getReturnType() {
        return CorrelationSymmetricMatrix.class;
    }
}
