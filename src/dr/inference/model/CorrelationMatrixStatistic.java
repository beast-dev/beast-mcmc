/*
 * MatrixInverseStatistic.java
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

package dr.inference.model;

import dr.math.matrixAlgebra.Matrix;
import dr.xml.*;

import static java.lang.Math.sqrt;

/**
 * @author Gabriel Hassler
 */

public class CorrelationMatrixStatistic extends Statistic.Abstract implements VariableListener {


    private static final String CORRELATION_MATRIX = "correlationMatrix";
    private final MatrixParameterInterface matrix;
    private final double[][] correlation;
    private boolean corrKnown = false;
    private final boolean invert;

    public CorrelationMatrixStatistic(MatrixParameterInterface matrix, Boolean invert) {
        this.matrix = matrix;
        this.invert = invert;
        correlation = new double[matrix.getRowDimension()][matrix.getColumnDimension()];
        matrix.addParameterListener(this);
    }


    @Override
    public int getDimension() {
        return matrix.getDimension();
    }

    @Override
    public double getStatisticValue(int dim) {
        if (!corrKnown) {

            double[][] varMat = matrix.getParameterAsMatrix();
            if (invert) {
                varMat = (new Matrix(varMat).inverse()).toComponents();
            }

            for (int i = 0; i < varMat.length; i++) {
                correlation[i][i] = 1;
                for (int j = i + 1; j < varMat.length; j++) {
                    double val = varMat[i][j] / sqrt(varMat[i][i] * varMat[j][j]);
                    correlation[i][j] = val;
                    correlation[j][i] = val;
                }
            }

            corrKnown = true;

        }
        int row = dim / correlation.length;
        int col = dim - row * correlation.length;

        return correlation[row][col];
    }

    @Override
    public String getDimensionName(int dim) {
        return getStatisticName() + "." + matrix.getDimensionName(dim);
    }


    @Override
    public void variableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
        corrKnown = false;
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        private static final String INVERT = "invert";

        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            MatrixParameterInterface matrix = (MatrixParameterInterface) xo.getChild(MatrixParameterInterface.class);
            if (matrix.getColumnDimension() != matrix.getRowDimension()) {
                throw new XMLParseException("Only square matrices can be converted to correlation matrices");
            }
            Boolean invert = false;
            if (xo.hasAttribute(INVERT)) {
                invert = xo.getBooleanAttribute(INVERT);
            }
            return new CorrelationMatrixStatistic(matrix, invert);
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{
                    new ElementRule(MatrixParameterInterface.class),
                    AttributeRule.newBooleanRule(INVERT, true)
            };
        }

        @Override
        public String getParserDescription() {
            return "This element returns a statistic that is the correlation matrix of the associated child statistic "
                    + "(potentially inverting the child statistic).";
        }

        @Override
        public Class getReturnType() {
            return CorrelationMatrixStatistic.class;
        }

        @Override
        public String getParserName() {
            return CORRELATION_MATRIX;
        }
    };
}
