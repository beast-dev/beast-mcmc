/*
 * MarkovRandomFieldMatrixParser.java
 *
 * Copyright (c) 2002-2012 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.inference.model.MarkovRandomFieldMatrix;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 *
 */
public class MarkovRandomFieldMatrixParser extends AbstractXMLObjectParser {

    public final static String MATRIX_PARAMETER = "markovRandomFieldMatrix";
    public static final String DIAGONAL = "diagonal";
    public static final String OFF_DIAGONAL = "offDiagonal";
    public static final String AS_CORRELATION = "asCorrelation";

    public String getParserName() {
        return MATRIX_PARAMETER;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        XMLObject cxo = xo.getChild(DIAGONAL);
        Parameter diagonalParameter = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(OFF_DIAGONAL);
        Parameter offDiagonalParameter = (Parameter) cxo.getChild(Parameter.class);

        boolean asCorrelation = xo.getAttribute(AS_CORRELATION, false);

        return new MarkovRandomFieldMatrix(diagonalParameter, offDiagonalParameter, asCorrelation);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A MRF matrix parameter constructed from its diagonals and first-order off diagonal.";
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(DIAGONAL,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(OFF_DIAGONAL,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            AttributeRule.newBooleanRule(AS_CORRELATION, true)
    };

    public Class getReturnType() {
        return MarkovRandomFieldMatrix.class;
    }
}
