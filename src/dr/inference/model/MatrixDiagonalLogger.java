/*
 * MatrixDiagonalLogger.java
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

public class MatrixDiagonalLogger extends Statistic.Abstract {
    private final MatrixParameterInterface matrix;

    public MatrixDiagonalLogger(MatrixParameterInterface matrix) {
        this.matrix = matrix;
    }


    @Override
    public int getDimension() {
        return matrix.getColumnDimension();
    }

    @Override
    public double getStatisticValue(int dim) {
        return matrix.getParameterValue(dim, dim);
    }

    @Override
    public String getDimensionName(int dim) {
        return getStatisticName() + "." + matrix.getDimensionName(dim);
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        private static final String MATRIX_DIAGONAL = "matrixDiagonals";


        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            MatrixParameterInterface matrix = (MatrixParameterInterface) xo.getChild(MatrixParameterInterface.class);
            if (matrix.getColumnDimension() != matrix.getRowDimension()) {
                throw new XMLParseException("Only square matrices can be converted to correlation matrices");
            }

            return new MatrixDiagonalLogger(matrix);
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{
                    new ElementRule(MatrixParameterInterface.class)
            };
        }

        @Override
        public String getParserDescription() {
            return "This element returns a statistic that is the diagonals of the associated matrix.";
        }

        @Override
        public Class getReturnType() {
            return MatrixDiagonalLogger.class;
        }

        @Override
        public String getParserName() {
            return MATRIX_DIAGONAL;
        }
    };
}
