/*
 * BlockDiagonalMatrixParameter.java
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

import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

//import java.util.StringTokenizer;

/**
 * @author Marc Suchard
 */
public class BlockDiagonalMatrixParameter extends MatrixParameter {

    public final static String BLOCK_DIAGONAL_MATRIX_PARAMETER = "blockDiagonalMatrixParameter";

    private final List<MatrixParameter> parameterList;

    public BlockDiagonalMatrixParameter(String name) {
        super(name);
        parameterList = new ArrayList<MatrixParameter>();
        rowOffset = new ArrayList<Integer>();
        colOffset = new ArrayList<Integer>();
    }

    private final List<Integer> rowOffset;
    private final List<Integer> colOffset;

//    public BlockDiagonalMatrixParameter(String name, Parameter[] parameters) {
//        super(name, parameters);
//        dimensionsEstablished = true;
//    }

    public double getParameterValue(int row, int col) {

        return 0; // TODO
    }

//    public double[][] getParameterAsMatrix() {
//        final int I = getRowDimension();
//        final int J = getColumnDimension();
//        double[][] parameterAsMatrix = new double[I][J];
//        for (int i = 0; i < I; i++) {
//            for (int j = 0; j < J; j++)
//                parameterAsMatrix[i][j] = getParameterValue(i, j);
//        }
//        return parameterAsMatrix;
//    }
//
//    public void setColumnDimension(int columnDimension) {
//        if (dimensionsEstablished) {
//            throw new IllegalArgumentException("Attempt to change dimensions after initialization");
//        }
//        this.columnDimension = columnDimension;
//        setupParameters();
//    }
//
//    public void setRowDimension(int rowDimension) {
//        if (dimensionsEstablished) {
//            throw new IllegalArgumentException("Attempt to change dimensions after initialization");
//        }
//        this.rowDimension = rowDimension;
//        setupParameters();
//    }
//
//    private void setupParameters() {
//        if (columnDimension > 0 && rowDimension > 0) {
//            dimensionsEstablished = true;
//
//            for (int i = 0; i < rowDimension; i++) {
//                Parameter row = new Parameter.Default(columnDimension, 0.0);
//                row.addBounds(new DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, columnDimension));
//                addParameter(row);
//            }
//        }
//    }

    public int getColumnDimension() {
        return columnDimension;
    }

    public int getRowDimension() {
        return rowDimension;
    }

//    public String toSymmetricString() {
//        StringBuilder sb = new StringBuilder("{");
//        int dim = getRowDimension();
//        int total = dim * (dim + 1) / 2;
//        for (int i = 0; i < dim; i++) {
//            for (int j = i; j < dim; j++) {
//                sb.append(String.format("%5.4e", getParameterValue(i, j)));
//                total--;
//                if (total > 0)
//                    sb.append(",");
//            }
//        }
//        sb.append("}");
//        return sb.toString();
//    }

//    public static BlockDiagonalMatrixParameter parseFromSymmetricString(String string) {
//        String clip = string.replace("{", "").replace("}", "").trim();
//        StringTokenizer st = new StringTokenizer(clip, ",");
//        int count = st.countTokens();
//        int dim = (-1 + (int) Math.sqrt(1 + 8 * count)) / 2;
//        Parameter[] parameter = new Parameter[dim];
//        for (int i = 0; i < dim; i++)
//            parameter[i] = new Parameter.Default(dim);
//        for (int i = 0; i < dim; i++) {
//            for (int j = i; j < dim; j++) {
//                double datum = new Double(st.nextToken());
//                parameter[i].setParameterValue(j, datum);
//                parameter[j].setParameterValue(i, datum);
//            }
//        }
//        return new BlockDiagonalMatrixParameter(null, parameter);
//    }
//
//    public static BlockDiagonalMatrixParameter parseFromSymmetricDoubleArray(Object[] data) {
//
//        int dim = (-1 + (int) Math.sqrt(1 + 8 * data.length)) / 2;
//        Parameter[] parameter = new Parameter[dim];
//        for (int i = 0; i < dim; i++)
//            parameter[i] = new Parameter.Default(dim);
//        int index = 0;
//        for (int i = 0; i < dim; i++) {
//            for (int j = i; j < dim; j++) {
//                double datum = (Double) data[index++];
//                parameter[i].setParameterValue(j, datum);
//                parameter[j].setParameterValue(i, datum);
//            }
//        }
//        return new BlockDiagonalMatrixParameter(null, parameter);
//    }

//    private boolean dimensionsEstablished = false;
    private int columnDimension = 0;
    private int rowDimension = 0;


    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return BLOCK_DIAGONAL_MATRIX_PARAMETER;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            final String name = xo.hasId() ? xo.getId() : null;

            BlockDiagonalMatrixParameter matrixParameter
                 = new BlockDiagonalMatrixParameter(name);

            for (int i = 0; i < xo.getChildCount(); i++) {
                MatrixParameter parameter = (MatrixParameter) xo.getChild(i);
                matrixParameter.addParameter(parameter); // TODO Double-check
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
                new ElementRule(MatrixParameter.class, 0, Integer.MAX_VALUE),
        };

        public Class getReturnType() {
            return BlockDiagonalMatrixParameter.class;
        }
    };
}
