/*
 * AdaptableSizeFastMatrixParameterParser.java
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


import dr.inference.model.AdaptableSizeFastMatrixParameter;
import dr.xml.*;

/**
 * Created by max on 4/6/16.
 */
public class AdaptableSizeFastMatrixParameterParser extends AbstractXMLObjectParser {


    private static final String ADAPTABLE_SIZE_FAST_MATRIX_PARAMETER="adaptableSizeFastMatrixParameter";
    private static final String ROWS="rows";
    private static final String MAX_ROW_SIZE="maxRowSize";
    private static final String MAX_COL_SIZE="maxColumnSize";
    private static final String COLUMNS="columns";
    private static final String TRANSPOSE="transpose";
    private static final String STARTING_VALUE = "startingValue";
    private static final String LOWER_TRIANGLE = "lowerTriangle";
    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        int MaxRowSize = xo.getAttribute(MAX_ROW_SIZE, 1);
        int MaxColumnSize = xo.getAttribute(MAX_COL_SIZE, 1);
        int rowDimension = xo.getAttribute(ROWS, 1);
        int columnDimension = xo.getAttribute(COLUMNS, 1);
        String name = xo.getId();
        double startingValue = 1;
        if(xo.hasAttribute(STARTING_VALUE))
            startingValue = xo.getDoubleAttribute(STARTING_VALUE);
        boolean lowerTriangle = false;
        if(xo.hasAttribute(LOWER_TRIANGLE))
            lowerTriangle = xo.getBooleanAttribute(LOWER_TRIANGLE);

      return new AdaptableSizeFastMatrixParameter(name, rowDimension, columnDimension, MaxRowSize, MaxColumnSize, startingValue, lowerTriangle);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
//            new ElementRule(Parameter.class, 0, Integer.MAX_VALUE),
            AttributeRule.newIntegerRule(COLUMNS, true),
            AttributeRule.newIntegerRule(ROWS, true),
            AttributeRule.newIntegerRule(MAX_ROW_SIZE, true),
            AttributeRule.newIntegerRule(MAX_COL_SIZE, true),
            AttributeRule.newDoubleRule(STARTING_VALUE, true),
            AttributeRule.newBooleanRule(LOWER_TRIANGLE, true),

    };



    @Override
    public String getParserDescription() {
        return "Returns a blockUpperTriangularMatrixParameter which is a compoundParameter which forces the last element to be of full length, the second to last element to be of full length-1, etc.";  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public Class getReturnType() {
        return AdaptableSizeFastMatrixParameter.class;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public String getParserName() {
        return ADAPTABLE_SIZE_FAST_MATRIX_PARAMETER;  //To change body of implemented methods use File | Settings | File Templates.
    }
}