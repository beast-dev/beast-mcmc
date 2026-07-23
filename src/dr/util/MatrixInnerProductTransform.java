/*
 * MatrixInnerProductTransform.java
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

package dr.util;

import dr.inference.model.MatrixParameterInterface;
import dr.xml.*;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

/**
 * @author Gabriel Hassler
 * @author Marc A. Suchard
 */

public class MatrixInnerProductTransform extends Transform.MatrixVariateTransform {


    private static final String MATRIX_INNER_PRODUCT = "matrixInnerProductTransform";


    public MatrixInnerProductTransform(int nRows, int nCols) {
        super(nRows * nCols, nRows, nCols);
    }

    @Override
    public double[] inverse(double[] values, int from, int to, double sum) {
        throw new RuntimeException(getTransformName() + " is not invertible");
    }

    @Override
    public double[] gradient(double[] values, int from, int to) {
        throw new RuntimeException("not yet implemented");
    }

    @Override
    public double[] gradientInverse(double[] values, int from, int to) {
        throw new RuntimeException("not yet implemented");
    }

    @Override
    public String getTransformName() {
        return MATRIX_INNER_PRODUCT;
    }


    @Override
    protected double[] transform(double[] values) {
        DenseMatrix64F X = DenseMatrix64F.wrap(rowDimension, columnDimension, values);
        DenseMatrix64F XXt = new DenseMatrix64F(rowDimension, rowDimension);
        CommonOps.multTransB(X, X, XXt);

        return XXt.getData();
    }

    @Override
    protected double[] inverse(double[] values) {
        throw new RuntimeException(getTransformName() + " is not invertible");
    }

    @Override
    protected double getLogJacobian(double[] values) {
        throw new RuntimeException(getTransformName() + " is not invertible");
    }

    @Override
    protected double[] getGradientLogJacobianInverse(double[] values) {
        throw new RuntimeException(getTransformName() + " is not invertible");
    }

    @Override
    public double[][] computeJacobianMatrixInverse(double[] values) {
        throw new RuntimeException(getTransformName() + " is not invertible");
    }

    @Override
    protected boolean isInInteriorDomain(double[] values) {
        return values.length == dim;
    }

    public static final AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {

        private static final String N_ROWS = "nRows";
        private static final String N_COLS = "nColumns";

        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            final int rowDimension;
            final int columnDimension;

            if (xo.getChildCount() > 0) {
                MatrixParameterInterface parameter = (MatrixParameterInterface)
                        xo.getChild(MatrixParameterInterface.class);

                rowDimension = parameter.getRowDimension();
                columnDimension = parameter.getColumnDimension();
            } else {
                rowDimension = xo.getIntegerAttribute(N_ROWS);
                columnDimension = xo.getIntegerAttribute(N_COLS);
            }

            return new MatrixInnerProductTransform(rowDimension, columnDimension);
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{
//                    new XORRule(
//                            new ElementRule(MatrixParameterInterface.class),
//                            new AndRule(
//                                    AttributeRule.newIntegerRule(N_ROWS),
//                                    AttributeRule.newIntegerRule(N_COLS))
//                    )
            };
        }

        @Override
        public String getParserDescription() {
            return "Takes the matrix X and transforms it to XXt";
        }

        @Override
        public Class getReturnType() {
            return MatrixInnerProductTransform.class;
        }

        @Override
        public String getParserName() {
            return MATRIX_INNER_PRODUCT;
        }
    };
}
