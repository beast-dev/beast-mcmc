/*
 * MaskedMatrixParameterParser.java
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

import dr.inference.model.MaskedParameter;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * @author Marc A. Suchard
 */
public class MaskedMatrixParameterParser extends AbstractXMLObjectParser {

    public static final String MASKED_MATRIX_PARAMETER = "maskedMatrixParameter";
    public static final String MASKING = "mask";
    public static final String COMPLEMENT = "complement";
    public static final String FROM = "from";
    public static final String TO = "to";
    public static final String EVERY = "every";

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        MatrixParameter matrix = (MatrixParameter) xo.getChild(MatrixParameter.class);


        System.err.println("colDim " + matrix.getColumnDimension());
        System.err.println("rowDim " + matrix.getRowDimension());
        System.err.println("parDim " + matrix.getParameterCount());


        Parameter mask;

        XMLObject cxo = xo.getChild(MASKING);
        if (cxo != null) {
            mask = (Parameter) cxo.getChild(Parameter.class);
        } else {
            int from = xo.getAttribute(FROM, 1) - 1;
            int to = xo.getAttribute(TO, matrix.getRowDimension()) - 1;
            int every = xo.getAttribute(EVERY, 1);

            if (from < 0) throw new XMLParseException("illegal 'from' attribute in " + MASKED_MATRIX_PARAMETER
                    + " element");
            if (to < 0 || to < from) throw new XMLParseException("illegal 'to' attribute in "
                    + MASKED_MATRIX_PARAMETER + " element");
            if (every <= 0) throw new XMLParseException("illegal 'every' attribute in " + MASKED_MATRIX_PARAMETER
                    + " element");

            mask = new Parameter.Default(matrix.getRowDimension(), 0.0);
            for (int i = from; i <= to; i += every) {
                mask.setParameterValue(i, 1.0);
            }
        }

        if (mask.getDimension() != matrix.getRowDimension())
            throw new XMLParseException("rowDim(" + matrix.getId() + ") != dim(" + mask.getId() + ")");

        boolean ones = !xo.getAttribute(COMPLEMENT, false);

        MaskedParameter[] maskedParameters = new MaskedParameter[matrix.getColumnDimension()];
        for (int col = 0; col < matrix.getColumnDimension(); ++col) {
            maskedParameters[col] = new MaskedParameter(matrix.getParameter(col));
            maskedParameters[col].addMask(mask, ones);
        }

        MatrixParameter maskedMatrix = new MatrixParameter(matrix.getId() + ".masked", maskedParameters);

//        for (int col = 0; col < matrix.getColumnDimension(); ++col) {
//            maskedMatrix.addParameter(matrix.getParameter(col));
//        }

        return maskedMatrix;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(MatrixParameter.class),
            new ElementRule(MASKING,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class)
                    }, true),
            AttributeRule.newBooleanRule(COMPLEMENT, true),
            AttributeRule.newIntegerRule(FROM, true),
            AttributeRule.newIntegerRule(TO, true),
            AttributeRule.newIntegerRule(EVERY, true),
    };

    public String getParserDescription() {
        return "A masked matrix parameter.";
    }

    public Class getReturnType() {
        return MatrixParameter.class;
    }

    public String getParserName() {
        return MASKED_MATRIX_PARAMETER;
    }
}
