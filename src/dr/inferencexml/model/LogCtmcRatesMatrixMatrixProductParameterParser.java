/*
 * LogCtmcRatesMatrixMatrixProductParameterParser.java
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

import dr.inference.model.LogCtmcRatesFromMatrixMatrixProductParameter;
import dr.inference.model.MatrixParameter;
import dr.xml.*;

/**
 * @author Marc A. Suchard
 */

public class LogCtmcRatesMatrixMatrixProductParameterParser extends AbstractXMLObjectParser {

    public static final String PRODUCT_PARAMETER = "ctmcRatesMatrixMatrixProduct";
    public static final String LHS = "lhs";
    public static final String RHS = "rhs";

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        MatrixParameter lhs = (MatrixParameter) xo.getChild(LHS).getChild(MatrixParameter.class);
        MatrixParameter rhs = (MatrixParameter) xo.getChild(RHS).getChild(MatrixParameter.class);

        if (lhs.getRowDimension() != rhs.getColumnDimension() ||
                lhs.getColumnDimension() != rhs.getRowDimension()) {
            throw new XMLParseException("Wrong matrix-matrix dimensions in " + xo.getId());
        }

        return new LogCtmcRatesFromMatrixMatrixProductParameter(lhs, rhs);

    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(LHS, new XMLSyntaxRule[]{
                    new ElementRule(MatrixParameter.class),
            }),
            new ElementRule(RHS, new XMLSyntaxRule[]{
                    new ElementRule(MatrixParameter.class),
            }),
    };

    public String getParserDescription() {
        return "A matrix-matrix product of parameters for CTMC rates.";
    }

    public Class getReturnType() {
        return LogCtmcRatesFromMatrixMatrixProductParameter.class;
    }

    public String getParserName() {
        return PRODUCT_PARAMETER;
    }
}
