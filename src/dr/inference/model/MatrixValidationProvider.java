/*
 * MatrixValidationProvider.java
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

package dr.inference.model;

import dr.xml.*;

public class MatrixValidationProvider implements CrossValidationProvider {

    private final MatrixParameter trueParameter;
    private final MatrixParameter inferredParameter;
    private final int[] relevantDimensions;
    private final String[] colNames;
    private final String sumName;

    MatrixValidationProvider(MatrixParameter trueParameter, MatrixParameter inferredParameter, String id) {
        this.trueParameter = trueParameter;
        this.inferredParameter = inferredParameter;

        int dimParameter = trueParameter.getDimension();

        this.relevantDimensions = new int[dimParameter];
        for (int i = 0; i < dimParameter; i++) {
            relevantDimensions[i] = i;
        }

        this.colNames = new String[dimParameter];

        for (int i = 0; i < dimParameter; i++) {

            int row = i / trueParameter.getRowDimension();
            int col = i - row * trueParameter.getRowDimension();

            colNames[i] = id + (row + 1) + (col + 1);

        }

        sumName = id + ".TotalSum";

    }

    @Override
    public double[] getTrueValues() {
        return trueParameter.getParameterValues();
    }

    @Override
    public double[] getInferredValues() {
        return inferredParameter.getParameterValues();
    }

    @Override
    public int[] getRelevantDimensions() {
        return relevantDimensions;
    }

    @Override
    public String getName(int dim) {
        return colNames[dim];
    }

    @Override
    public String getNameSum(int dim) {
        return sumName;
    }

    //TODO: Merge with TraitValidationProvider parser ?
    public static dr.xml.XMLObjectParser PARSER = new dr.xml.AbstractXMLObjectParser() {

        final static String PARSER_NAME = "matrixValidation";
        final static String TRUE_PARAMETER = "trueParameter";
        final static String INFERRED_PARAMETER = "inferredParameter";
        final static String LOG_SUM = "logSum";

        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            MatrixParameter trueParameter = (MatrixParameter) xo.getElementFirstChild(TRUE_PARAMETER);
            MatrixParameter inferredParameter = (MatrixParameter) xo.getElementFirstChild(INFERRED_PARAMETER);


            String id = PARSER_NAME;
            if (xo.hasId()) {
                id = xo.getId();
            }

            if (trueParameter.getRowDimension() != inferredParameter.getRowDimension() ||
                    trueParameter.getColumnDimension() != inferredParameter.getColumnDimension()) {

                throw new XMLParseException("The matrix parameters contained in " + TRUE_PARAMETER + " and " + INFERRED_PARAMETER +
                        " must have the same dimensions.");
            }

            MatrixValidationProvider provider = new MatrixValidationProvider(trueParameter, inferredParameter, id);

            boolean logSum = xo.getAttribute(LOG_SUM, false);

            if (logSum) return new CrossValidatorSum(provider, ValidationType.SQUARED_ERROR);
            return new CrossValidator(provider, ValidationType.SQUARED_ERROR);

        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{
                    AttributeRule.newBooleanRule(LOG_SUM, true),
                    new ElementRule(TRUE_PARAMETER, new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class)
                    }),
                    new ElementRule(INFERRED_PARAMETER, new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class)
                    })
            };
        }

        @Override
        public String getParserDescription() {
            return null;
        }

        @Override
        public Class getReturnType() {
            return CrossValidator.class;
        }

        @Override
        public String getParserName() {
            return PARSER_NAME;
        }


    };
}
