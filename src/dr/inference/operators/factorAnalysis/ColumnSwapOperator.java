/*
 * ColumnSwapOperator.java
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

package dr.inference.operators.factorAnalysis;

import dr.inference.model.MatrixParameterInterface;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;
import dr.xml.*;

public class ColumnSwapOperator extends SimpleMCMCOperator {

    private final MatrixParameterInterface matrix;
    private static final String COLUMN_SWAP = "columnSwapOperator";

    public ColumnSwapOperator(MatrixParameterInterface matrix, double weight) {
        setWeight(weight);
        this.matrix = matrix;
    }


    @Override
    public String getOperatorName() {
        return COLUMN_SWAP;
    }

    @Override
    public double doOperation() {
        int n = matrix.getRowDimension();
        int p = matrix.getColumnDimension();

        int i1 = MathUtils.nextInt(p);
        int i2 = MathUtils.nextInt(p);
        while (i2 == i1) {
            i2 = MathUtils.nextInt(p);
        }

        for (int i = 0; i < n; i++) {
            double v1 = matrix.getParameterValue(i, i1);
            double v2 = matrix.getParameterValue(i, i2);
            matrix.setParameterValueQuietly(i, i1, v2);
            matrix.setParameterValueQuietly(i, i2, v1);
        }
        matrix.fireParameterChangedEvent();
        return 0;
    }

    private static final String WEIGHT = "weight";

    public static AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {
        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            double weight = xo.getDoubleAttribute(WEIGHT);
            MatrixParameterInterface matrix = (MatrixParameterInterface) xo.getChild(MatrixParameterInterface.class);

            return new ColumnSwapOperator(matrix, weight);
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{
                    new ElementRule(MatrixParameterInterface.class),
                    AttributeRule.newDoubleRule(WEIGHT)
            };
        }

        @Override
        public String getParserDescription() {
            return "swaps random columns of a matrix";
        }

        @Override
        public Class getReturnType() {
            return ColumnSwapOperator.class;
        }

        @Override
        public String getParserName() {
            return "columnSwapOperator";
        }
    };
}
