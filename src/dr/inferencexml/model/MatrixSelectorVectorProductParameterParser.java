/*
 * MatrixSelectorVectorProductParameterParser.java
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

package dr.inferencexml.model;

import dr.inference.model.MatrixParameter;
import dr.inference.model.MatrixSelectorVectorProductParameter;
import dr.inference.model.MatrixVectorProductParameter;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * @author Marc A. Suchard
 */

public class MatrixSelectorVectorProductParameterParser extends AbstractXMLObjectParser {

    public static final String SELECTOR_PRODUCT_PARAMETER = "matrixSelectorVectorProductParameter";
    public static final String MATRIX = "matrix";
    public static final String VECTOR = "vector";
    private static final String SELECTOR = "selector";

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        MatrixParameter matrix = (MatrixParameter) xo.getChild(MATRIX).getChild(MatrixParameter.class);
        Parameter vector = (Parameter) xo.getChild(VECTOR).getChild(Parameter.class);
        Parameter selector = (Parameter) xo.getChild(SELECTOR).getChild(Parameter.class);

        if (matrix.getColumnDimension() != selector.getDimension()) {
            throw new XMLParseException("Wrong matrix-selector dimensions in " + xo.getId());
        }

        return new MatrixSelectorVectorProductParameter(matrix, vector, selector);

    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(MATRIX, new XMLSyntaxRule[]{
                    new ElementRule(MatrixParameter.class),
            }),
            new ElementRule(VECTOR, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class),
            }),
    };

    public String getParserDescription() {
        return "A matrix-selector-vector product of parameters.";
    }

    public Class getReturnType() {
        return MatrixSelectorVectorProductParameter.class;
    }

    public String getParserName() {
        return SELECTOR_PRODUCT_PARAMETER;
    }
}
