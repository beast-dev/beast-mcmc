/*
 * MatrixParameter.java
 *
 * Copyright (c) 2002-2012 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.xml.*;

import java.util.StringTokenizer;

/**
 * @author Marc Suchard
 */
public class MatrixParameter extends CompoundParameter {

    public final static String MATRIX_PARAMETER = "matrixParameter";

    public MatrixParameter(String name) {
        super(name);
    }

    public MatrixParameter(String name, Parameter[] parameters) {
        super(name, parameters);
        dimensionsEstablished = true;
    }

    public double getParameterValue(int row, int col) {
        return getParameter(col).getParameterValue(row);
    }

    public double[][] getParameterAsMatrix() {
        final int I = getRowDimension();
        final int J = getColumnDimension();
        double[][] parameterAsMatrix = new double[I][J];
        for (int i = 0; i < I; i++) {
            for (int j = 0; j < J; j++)
                parameterAsMatrix[i][j] = getParameterValue(i, j);
        }
        return parameterAsMatrix;
    }

    public void setColumnDimension(int columnDimension) {
        if (dimensionsEstablished) {
            throw new IllegalArgumentException("Attempt to change dimensions after initialization");
        }
        this.columnDimension = columnDimension;
        setupParameters();
    }

    public void setRowDimension(int rowDimension) {
        if (dimensionsEstablished) {
            throw new IllegalArgumentException("Attempt to change dimensions after initialization");
        }
        this.rowDimension = rowDimension;
        setupParameters();
    }

    private void setupParameters() {
        if (columnDimension > 0 && rowDimension > 0) {
            dimensionsEstablished = true;

            for (int i = 0; i < rowDimension; i++) {
                Parameter row = new Parameter.Default(columnDimension, 0.0);
                row.addBounds(new DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, columnDimension));
                addParameter(row);
            }
        }
    }

    public int getColumnDimension() {
        return getParameterCount();
    }

    public int getRowDimension() {
        return getParameter(0).getDimension();
    }

    public String toSymmetricString() {
        StringBuilder sb = new StringBuilder("{");
        int dim = getRowDimension();
        int total = dim * (dim + 1) / 2;
        for (int i = 0; i < dim; i++) {
            for (int j = i; j < dim; j++) {
                sb.append(String.format("%5.4e", getParameterValue(i, j)));
                total--;
                if (total > 0)
                    sb.append(",");
            }
        }
        sb.append("}");
        return sb.toString();
    }

    public static MatrixParameter parseFromSymmetricString(String string) {
        String clip = string.replace("{", "").replace("}", "").trim();
        StringTokenizer st = new StringTokenizer(clip, ",");
        int count = st.countTokens();
        int dim = (-1 + (int) Math.sqrt(1 + 8 * count)) / 2;
        Parameter[] parameter = new Parameter[dim];
        for (int i = 0; i < dim; i++)
            parameter[i] = new Parameter.Default(dim);
        for (int i = 0; i < dim; i++) {
            for (int j = i; j < dim; j++) {
                double datum = new Double(st.nextToken());
                parameter[i].setParameterValue(j, datum);
                parameter[j].setParameterValue(i, datum);
            }
        }
        return new MatrixParameter(null, parameter);
    }

    public static MatrixParameter parseFromSymmetricDoubleArray(Object[] data) {

        int dim = (-1 + (int) Math.sqrt(1 + 8 * data.length)) / 2;
        Parameter[] parameter = new Parameter[dim];
        for (int i = 0; i < dim; i++)
            parameter[i] = new Parameter.Default(dim);
        int index = 0;
        for (int i = 0; i < dim; i++) {
            for (int j = i; j < dim; j++) {
                double datum = (Double) data[index++];
                parameter[i].setParameterValue(j, datum);
                parameter[j].setParameterValue(i, datum);
            }
        }
        return new MatrixParameter(null, parameter);
    }

    private boolean dimensionsEstablished = false;
    private int columnDimension = 0;
    private int rowDimension = 0;

    // **************************************************************
    // XMLElement IMPLEMENTATION
    // **************************************************************

//    public Element createElement(Document d) {
//        throw new RuntimeException("Not implemented yet!");
//    }

    private static final String ROW_DIMENSION = "rows";
    private static final String COLUMN_DIMENSION = "columns";
    private static final String TRANSPOSE = "transpose";

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return MATRIX_PARAMETER;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            final String name = xo.hasId() ? xo.getId() : null;
            boolean transposed = xo.getAttribute(TRANSPOSE, false);

            MatrixParameter matrixParameter;
            if (!transposed) {
                matrixParameter = new MatrixParameter(name);
            } else {
                matrixParameter = new TransposedMatrixParameter(name);
            }

            if (xo.hasAttribute(ROW_DIMENSION)) {
                int rowDimension = xo.getIntegerAttribute(ROW_DIMENSION);
                matrixParameter.setRowDimension(rowDimension);
            }

            if (xo.hasAttribute(COLUMN_DIMENSION)) {
                int columnDimension = xo.getIntegerAttribute(COLUMN_DIMENSION);
                matrixParameter.setColumnDimension(columnDimension);
            }

            int dim = 0;
            for (int i = 0; i < xo.getChildCount(); i++) {
                Parameter parameter = (Parameter) xo.getChild(i);
                matrixParameter.addParameter(parameter);
                if (i == 0)
                    dim = parameter.getDimension();
                else if (dim != parameter.getDimension())
                    throw new XMLParseException("All parameters must have the same dimension to construct a rectangular matrix");
            }

            return matrixParameter;
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "A matrix parameter constructed from its component parameters.";
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                new ElementRule(Parameter.class, 0, Integer.MAX_VALUE),
                AttributeRule.newIntegerRule(ROW_DIMENSION, true),
                AttributeRule.newIntegerRule(COLUMN_DIMENSION, true)
        };

        public Class getReturnType() {
            return MatrixParameter.class;
        }
    };
}
