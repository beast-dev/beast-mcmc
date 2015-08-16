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

/**
 * @author Marc A. Suchard
 */
public class MatrixInverseStatistic extends Statistic.Abstract implements VariableListener {

    public static final String INVERSE_STATISTIC = "matrixInverse";

    public MatrixInverseStatistic(MatrixParameter matrix) {
        this.matrix = matrix;
        matrix.addParameterListener(this);
    }

    public int getDimension() {
        return matrix.getDimension();
    }

    public double getStatisticValue(int dim) {

        if (!inverseKnown) {
            inverse = (new Matrix(matrix.getParameterAsMatrix()).inverse().toComponents());
            inverseKnown = true;
        }

        int x = dim / matrix.getColumnDimension();
        int y = dim - x * matrix.getColumnDimension();

        return inverse[x][y];
    }

    public String getDimensionName(int dim) {        
        return getStatisticName() + "." + matrix.getDimensionName(dim);
    }

    public void variableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        inverseKnown = false;
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return INVERSE_STATISTIC;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            MatrixParameter matrix = (MatrixParameter) xo.getChild(MatrixParameter.class);

            if (matrix.getColumnDimension() != matrix.getRowDimension())
                throw new XMLParseException("Only square matrices can be inverted");

            return new MatrixInverseStatistic(matrix);

        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element returns a statistic that is the matrix inverse of the child statistic.";
        }

        public Class getReturnType() {
            return MatrixInverseStatistic.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(MatrixParameter.class)
        };
    };

    private boolean inverseKnown = false;
    private double[][] inverse;
    private MatrixParameter matrix;
}
