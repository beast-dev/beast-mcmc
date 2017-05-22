///*
// * NewMatrixParameter.java
// *
// * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
// *
// * This file is part of BEAST.
// * See the NOTICE file distributed with this work for additional
// * information regarding copyright ownership and licensing.
// *
// * BEAST is free software; you can redistribute it and/or modify
// * it under the terms of the GNU Lesser General Public License as
// * published by the Free Software Foundation; either version 2
// * of the License, or (at your option) any later version.
// *
// *  BEAST is distributed in the hope that it will be useful,
// *  but WITHOUT ANY WARRANTY; without even the implied warranty of
// *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// *  GNU Lesser General Public License for more details.
// *
// * You should have received a copy of the GNU Lesser General Public
// * License along with BEAST; if not, write to the
// * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
// * Boston, MA  02110-1301  USA
// */
//
//package dr.inference.model;
//
//import dr.xml.*;
//
//import java.util.StringTokenizer;
////import org.w3c.dom.Document;
////import org.w3c.dom.Element;
////import dr.xml.*;
//
///**
// * @author Marc Suchard
// */
//public class NewMatrixParameter extends Parameter.Default {
//
//    public final static String MATRIX_PARAMETER = "matrixParameter";
//
//    public NewMatrixParameter(String name) {
//        super(name);
//    }
//
//    public NewMatrixParameter(String name, double[] parameter, int numRows, int numCols) {
//        super(parameter);
//        this.numRows = numRows;
//        this.numCols = numCols;
//    }
//
//    private int numRows = 0;
//    private int numCols = 0;
//
//    public double getParameterValue(int row, int col) {
//        return getParameterValue(row * numCols + col); // Stores in row-major
//    }
//
//    public double[][] getParameterAsMatrix() {
////        final int I = getRowDimension();
////        final int J = getColumnDimension();
////        double[][] parameterAsMatrix = new double[I][J];
////        for (int i = 0; i < I; i++) {
////            for (int j = 0; j < J; j++)
////                parameterAsMatrix[i][j] = getParameterValue(i, j);
////        }
////        return parameterAsMatrix;
//        double[][] parameterAsMatris = new double[numRows][numCols];
//        for (int i = 0; i < numRows; i++) {
//            for (int j = 0; j < numCols; j++) {
//                parameterAsMatris[i][j] = getParameterValue(i * numCols + j);
//            }
//        }
//        return parameterAsMatris;
//    }
//
//    public int getColumnDimension() {
//        return numCols;
//    }
//
//    public int getRowDimension() {
//        return numRows;
//    }
//
//    public String toSymmetricString() {
//        StringBuffer sb = new StringBuffer("{");
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
//
//    public static NewMatrixParameter parseFromSymmetricString(String string) {
//        String clip = string.replace("{", "").replace("}", "").trim();
//        StringTokenizer st = new StringTokenizer(clip, ",");
//        int count = st.countTokens();
//        int dim = (-1 + (int) Math.sqrt(1 + 8 * count)) / 2;
////        Parameter[] parameter = new Parameter[dim];
////        for (int i = 0; i < dim; i++)
////            parameter[i] = new Parameter.Default(dim);
////        for (int i = 0; i < dim; i++) {
////            for (int j = i; j < dim; j++) {
////                double datum = new Double(st.nextToken());
////                parameter[i].setParameterValue(j, datum);
////                parameter[j].setParameterValue(i, datum);
////            }
////        }
//        double[] data = new double[dim * dim];
//        for (int i = 0; i < dim; i++) {
//            for (int j = i; j < dim; j++) {
//                double datum = new Double(st.nextToken());
//                data[i * dim + j] = datum;
//                data[j * dim + i] = datum;
//            }
//        }
//        return new NewMatrixParameter(null, data, dim, dim);
//    }
//
//    public static NewMatrixParameter parseFromSymmetricDoubleArray(Object[] inData) {
//
//        int dim = (-1 + (int) Math.sqrt(1 + 8 * inData.length)) / 2;
////        Parameter[] parameter = new Parameter[dim];
////        for (int i = 0; i < dim; i++)
////            parameter[i] = new Parameter.Default(dim);
//        int index = 0;
////        for (int i = 0; i < dim; i++) {
////            for (int j = i; j < dim; j++) {
////                double datum = (Double) data[index++];
////                parameter[i].setParameterValue(j, datum);
////                parameter[j].setParameterValue(i, datum);
////            }
////        }
////        return new MatrixParameter(null, parameter);
//        double[] data = new double[dim * dim];
//        for (int i = 0; i < dim; i++) {
//            for (int j = i; j < dim; j++) {
//                double datum = (Double) inData[index++];
//                data[i * dim + j] = datum;
//                data[j * dim + i] = datum;
//            }
//        }
//        return new NewMatrixParameter(null, data, dim, dim);
//    }
//
//    // **************************************************************
//    // XMLElement IMPLEMENTATION
//    // **************************************************************
//
////    public Element createElement(Document d) {
////        throw new RuntimeException("Not implemented yet!");
////    }
//
//    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {
//
//        public String getParserName() {
//            return MATRIX_PARAMETER;
//        }
//
//        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
//
////            MatrixParameter matrixParameter = new MatrixParameter(MATRIX_PARAMETER);
//
//            int numRows = xo.getChildCount();
//            int numCols = 0;
//
//            double[] values = null;
//
//            for (int i = 0; i < numRows; i++) {
//                Parameter parameter = (Parameter) xo.getChild(i);
//                if (values == null) {
//                    numCols = parameter.getDimension();
//                    values = new double[numCols * numRows];
//                } else {
//                    if (numCols != parameter.getDimension()) {
//                        throw new XMLParseException(
//                                "All parameters must have the same dimension to construct a rectangular matrix");
//                    }
//                }
//                double[] newValues = parameter.getParameterValues();
//                System.arraycopy(newValues, 0, values, i * numCols, numCols);
//            }
//
//            String name = (xo.hasId() ? xo.getId() : MATRIX_PARAMETER);
//
//            return new NewMatrixParameter(name, values, numRows, numCols);
//        }
//
//        //************************************************************************
//        // AbstractXMLObjectParser implementation
//        //************************************************************************
//
//        public String getParserDescription() {
//            return "A matrix parameter constructed from its component parameters.";
//        }
//
//        public XMLSyntaxRule[] getSyntaxRules() {
//            return rules;
//        }
//
//        private final XMLSyntaxRule[] rules = {
//                new ElementRule(Parameter.class, 1, Integer.MAX_VALUE),
//        };
//
//        public Class getReturnType() {
//            return NewMatrixParameter.class;
//        }
//    };
//
//
//}