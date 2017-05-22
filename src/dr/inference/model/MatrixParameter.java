/*
 * MatrixParameter.java
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

import java.util.StringTokenizer;

/**
 * @author Marc Suchard
 * @author Max Tolkoff
 */
public class MatrixParameter extends CompoundParameter implements MatrixParameterInterface {

    public final static String MATRIX_PARAMETER = "matrixParameter";

    public MatrixParameter(String name) {
        super(name);
    }

    public MatrixParameter(String name, Parameter[] parameters) {
        super(name, parameters);
        rowDimension=parameters[0].getDimension();
        columnDimension=parameters.length;
        dimensionsEstablished = true;
    }

    public MatrixParameter(String name, int row, int column){
        super(name);
        setDimensions(row, column);
    }

    public MatrixParameter(String name, int row, int column, double a){
        super(name);
        setDimensions(row, column, a);
    }

    public void setParameterValue(int row, int column, double a) {
        getParameter(column).setParameterValue(row, a);
    }

    public void setParameterValueQuietly(int row, int column, double a){
        getParameter(column).setParameterValueQuietly(row, a);
    }

    public void setParameterValueNotifyChangedAll(int row, int column, double val){
        getParameter(column).setParameterValueNotifyChangedAll(row, val);
    }

    public static MatrixParameter recast(String name, CompoundParameter compoundParameter) {
        final int count = compoundParameter.getParameterCount();
        Parameter[] parameters = new Parameter[count];
        for (int i = 0; i < count; ++i) {
            parameters[i] = compoundParameter.getParameter(i);
        }
        return new MatrixParameter(name, parameters);
    }

    public double getParameterValue(int row, int col) {
        return getParameter(col).getParameterValue(row);
    }

    public double[] getRowValues(int row){
        int colDim=getColumnDimension();
        double[] rowValues=new double[colDim];
        for (int i = 0; i <colDim ; i++) {
            rowValues[i]=getParameterValue(row, i);
        }
        return rowValues;
    }

    public double[] getColumnValues(int col) {
        return this.getParameter(col).getParameterValues();
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


    //TODO rewrite so that it doesn't destroy existing parameters
    public void setDimensions(int rowDim, int colDim){
        setDimensions(rowDim, colDim, 0.0);
    }

    public void setDimensions(int rowDim, int colDim, double a){
        rowDimension=rowDim;
        columnDimension=colDim;
        for (int i = 0; i < colDim; i++) {
            Parameter column = new Parameter.Default(rowDim, a);
            column.addBounds(new DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, rowDim));
            addParameter(column);
        }
    }

    public int getColumnDimension() {
        return getParameterCount();
    }

    public int getRowDimension() {
        return getParameter(0).getDimension();
    }

    @Override
    public int getUniqueParameterCount() {
        return getParameterCount();
    }

    @Override
    public Parameter getUniqueParameter(int index) {
        return super.getParameter(index);
    }

    @Override
    public void copyParameterValues(double[] destination, int offset) {
        final int length = getDimension();
        for (int i = 0; i < length; ++i) {
            destination[offset + i] = getParameterValue(i);
        }
    }

    @Override
    public void setAllParameterValuesQuietly(double[] values, int offset) {
        for (int i = 0; i < getDimension(); ++i) {
            setParameterValueQuietly(i, values[offset + i]);
        }
    }

    public static String toSymmetricString(MatrixParameterInterface mat) {
        StringBuilder sb = new StringBuilder("{");
        int dim = mat.getRowDimension();
        int total = dim * (dim + 1) / 2;
        for (int i = 0; i < dim; i++) {
            for (int j = i; j < dim; j++) {
                sb.append(String.format("%5.4e", mat.getParameterValue(i, j)));
                total--;
                if (total > 0)
                    sb.append(",");
            }
        }
        sb.append("}");
        return sb.toString();
    }

    public String toSymmetricString() {
        return toSymmetricString(this);
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

    public void rowMultiply(double a, int row){
        rowMultiplyQuietly(a, row);
        fireParameterChangedEvent();
    }

    public void columnMultiply(double a, int col){
        columnMultiplyQuietly(a,col);
        fireParameterChangedEvent();
    }





    public TransposedMatrixParameter transpose(){
        return TransposedMatrixParameter.recast(null, this);
    }

    // **************************************************************
    // Matrix Operations
    // **************************************************************

    public MatrixParameter add(MatrixParameter Right) {
        if (Right.getRowDimension() != getRowDimension() || getColumnDimension() != Right.getColumnDimension()){
            throw new RuntimeException("You cannot add a " + getRowDimension() +" by " + getColumnDimension() + " matrix to a " + Right.getRowDimension() + " by " + Right.getColumnDimension() + " matrix.");
        }
        MatrixParameter answer=new MatrixParameter(null);
        answer.setDimensions(getRowDimension(), getColumnDimension());
        for (int i = 0; i <getRowDimension() ; i++) {
            for (int j = 0; j <getColumnDimension() ; j++) {
                answer.setParameterValueQuietly(i, j, getParameterValue(i, j) + Right.getParameterValue(i, j));
            }

        }
        return answer;
    }

    public MatrixParameter addInPlace(MatrixParameter Right, MatrixParameter answer) {
        if (Right.getRowDimension() != getRowDimension() || getColumnDimension() != Right.getColumnDimension()){
            throw new RuntimeException("You cannot add a " + getRowDimension() +" by " + getColumnDimension() + " matrix to a " + Right.getRowDimension() + " by " + Right.getColumnDimension() + " matrix.");
        }
//        MatrixParameter answer=new MatrixParameter(null);
//        answer.setDimensions(getRowDimension(), getColumnDimension());
        for (int i = 0; i <getRowDimension() ; i++) {
            for (int j = 0; j <getColumnDimension() ; j++) {
                answer.setParameterValueQuietly(i, j, getParameterValue(i, j) + Right.getParameterValue(i, j));
            }

        }
        return answer;
    }

    public MatrixParameter subtract(MatrixParameter Right) {
        if (Right.getRowDimension() != getRowDimension() || getColumnDimension() != Right.getColumnDimension()){
            throw new RuntimeException("You cannot subtract a " + getRowDimension() +" by " + getColumnDimension() + " matrix to a " + Right.getRowDimension() + " by " + Right.getColumnDimension() + " matrix.");
        }
        MatrixParameter answer=new MatrixParameter(null);
        answer.setDimensions(getRowDimension(), getColumnDimension());
        for (int i = 0; i <getRowDimension() ; i++) {
            for (int j = 0; j <getColumnDimension() ; j++) {
                answer.setParameterValueQuietly(i, j, getParameterValue(i, j) - Right.getParameterValue(i, j));
            }

        }
        return answer;
    }

    public MatrixParameter subtractInPlace(MatrixParameter Right, MatrixParameter answer) {
        if (Right.getRowDimension() != getRowDimension() || getColumnDimension() != Right.getColumnDimension()){
            throw new RuntimeException("You cannot subtract a " + getRowDimension() +" by " + getColumnDimension() + " matrix from a " + Right.getRowDimension() + " by " + Right.getColumnDimension() + " matrix.");
        }
//        MatrixParameter answer=new MatrixParameter(null);
//        answer.setDimensions(getRowDimension(), getColumnDimension());
        for (int i = 0; i <getRowDimension() ; i++) {
            for (int j = 0; j <getColumnDimension() ; j++) {
                answer.setParameterValueQuietly(i,j, getParameterValue(i,j)-Right.getParameterValue(i,j));
            }

        }
        return answer;
    }

    public MatrixParameter transposeThenProduct(MatrixParameter Right){
        if(this.getRowDimension()!=Right.getRowDimension()){
            throw new RuntimeException("Incompatible Dimensions: " + Right.getRowDimension() + " does not equal " + this.getRowDimension() +".\n");
        }
        MatrixParameter answer=new MatrixParameter(null);
        answer.setDimensions(this.getColumnDimension(), Right.getColumnDimension());

        int p = this.getRowDimension();
        int n = this.getColumnDimension();
        int m = Right.getColumnDimension();
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                double sum = 0;
                for (int k = 0; k < p; k++)
                    sum += getParameterValue(k,i) * Right.getParameterValue(k,j);
                answer.setParameterValueQuietly(i, j, sum);
            }
        }

        return answer;
    }

    public MatrixParameter transposeThenProductInPlace(MatrixParameter Right, MatrixParameter answer){
        if(this.getRowDimension()!=Right.getRowDimension()){
            throw new RuntimeException("Incompatible Dimensions: " + Right.getRowDimension() + " does not equal " + this.getRowDimension() +".\n");
        }
//        MatrixParameter answer=new MatrixParameter(null);
//        answer.setDimensions(this.getColumnDimension(), Right.getColumnDimension());

        int p = this.getRowDimension();
        int n = this.getColumnDimension();
        int m = Right.getColumnDimension();
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                double sum = 0;
                for (int k = 0; k < p; k++)
                    sum += getParameterValue(k,i) * Right.getParameterValue(k,j);
                answer.setParameterValueQuietly(i,j, sum);
            }
        }

        return answer;
    }


    public MatrixParameter product(MatrixParameter Right){
        if(this.getColumnDimension()!=Right.getRowDimension()){
            throw new RuntimeException("Incompatible Dimensions: " + Right.getRowDimension() + " does not equal " + this.getColumnDimension() +".\n");
        }
        MatrixParameter answer=new MatrixParameter(null);
        answer.setDimensions(this.getRowDimension(), Right.getColumnDimension());

        int p = this.getColumnDimension();
        int n = this.getRowDimension();
        int m = Right.getColumnDimension();
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                double sum = 0;
                for (int k = 0; k < p; k++)
                    sum += getParameterValue(i,k) * Right.getParameterValue(k,j);
                answer.setParameterValueQuietly(i, j, sum);
            }
        }

        return answer;
    }

    public MatrixParameter productInPlace(MatrixParameter Right, MatrixParameter answer){
        if(this.getColumnDimension()!=Right.getRowDimension()){
            throw new RuntimeException("Incompatible Dimensions: " + Right.getRowDimension() + " does not equal " + this.getColumnDimension() +".\n");
        }
//        MatrixParameter answer=new MatrixParameter(null);
//        answer.setDimensions(this.getRowDimension(), Right.getColumnDimension());

        int p = this.getColumnDimension();
        int n = this.getRowDimension();
        int m = Right.getColumnDimension();
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                double sum = 0;
                for (int k = 0; k < p; k++)
                    sum += getParameterValue(i,k) * Right.getParameterValue(k,j);
                answer.setParameterValueQuietly(i,j, sum);
            }
        }

        return answer;
    }

    public MatrixParameter productWithTransposed(MatrixParameter Right){
        if(this.getColumnDimension()!=Right.getColumnDimension()){
            throw new RuntimeException("Incompatible Dimensions: " + Right.getColumnDimension() + " does not equal " + this.getColumnDimension() +".\n");
        }
        MatrixParameter answer=new MatrixParameter(null);
        answer.setDimensions(this.getRowDimension(), Right.getRowDimension());

        int p = this.getColumnDimension();
        int n = this.getRowDimension();
        int m = Right.getRowDimension();
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                double sum = 0;
                for (int k = 0; k < p; k++)
                    sum += getParameterValue(i,k) * Right.getParameterValue(j,k);
                answer.setParameterValueQuietly(i,j, sum);
            }
        }

        return answer;
    }

    public MatrixParameter productWithTransposedInPlace(MatrixParameter Right, MatrixParameter answer){
        if(this.getColumnDimension()!=Right.getColumnDimension()){
            throw new RuntimeException("Incompatible Dimensions: " + Right.getColumnDimension() + " does not equal " + this.getColumnDimension() +".\n");
        }
//        MatrixParameter answer=new MatrixParameter(null);
//        answer.setDimensions(this.getRowDimension(), Right.getRowDimension());

        int p = this.getColumnDimension();
        int n = this.getRowDimension();
        int m = Right.getRowDimension();
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                double sum = 0;
                for (int k = 0; k < p; k++)
                    sum += getParameterValue(i,k) * Right.getParameterValue(j,k);
                answer.setParameterValueQuietly(i,j, sum);
            }
        }

        return answer;
    }

    public void product(double a){
        for (int i = 0; i <this.getRowDimension() ; i++) {
            for (int j = 0; j < this.getColumnDimension(); j++) {
                this.setParameterValueQuietly(i, j, a * this.getParameterValue(i, j));
            }
        }
        this.fireParameterChangedEvent();
    }

    public MatrixParameter productInPlace(double a, MatrixParameter answer){
//        MatrixParameter answer=new MatrixParameter(null);
//        answer.setDimensions(this.getRowDimension(), this.getColumnDimension());
        for (int i = 0; i <this.getRowDimension() ; i++) {
            for (int j = 0; j <this.getColumnDimension() ; j++) {
                answer.setParameterValueQuietly(i,j, a*this.getParameterValue(i,j));
            }

        }
        return answer;
    }

    public void rowMultiplyQuietly(double a, int row){
        for (int i = 0; i <getColumnDimension() ; i++) {
            if(getParameter(i).getDimension()<row){
                getParameter(i).setParameterValueQuietly(row, a*getParameterValue(row, i));
            }
        }
    }

    public void columnMultiplyQuietly(double a, int col){
        Parameter i=getParameter(col);
        for (int j = 0; j < i.getDimension() ; j++) {
            i.setParameterValueQuietly(j, a*i.getParameterValue(j));
        }
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

    public static final String ROW_DIMENSION = "rows";
    public static final String COLUMN_DIMENSION = "columns";
    public static final String TRANSPOSE = "transpose";
    public static final String AS_COMPOUND = "asCompoundParameter";
    public static final String BEHAVIOR = "test";

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return MATRIX_PARAMETER;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            final String name = xo.hasId() ? xo.getId() : null;
            boolean transposed = xo.getAttribute(TRANSPOSE, false);
            boolean compound = xo.getAttribute(AS_COMPOUND, false);

            MatrixParameter matrixParameter;

            if (compound) {
                CompoundParameter parameter = (CompoundParameter) xo.getChild(0);
                if (transposed) {
                    matrixParameter = TransposedMatrixParameter.recast(name, parameter);
                } else {
                    matrixParameter = MatrixParameter.recast(name, parameter);
                }
                return matrixParameter;
            }

            if (!transposed) {
                matrixParameter = new MatrixParameter(name);
            } else {
                matrixParameter = new TransposedMatrixParameter(name);
            }

            if (xo.getAttribute(BEHAVIOR, false) && xo.hasAttribute(ROW_DIMENSION) && xo.hasAttribute(COLUMN_DIMENSION)) {
                int rowDim = xo.getIntegerAttribute(ROW_DIMENSION);
                int colDim = xo.getIntegerAttribute(COLUMN_DIMENSION);
                matrixParameter.setDimensions(rowDim, colDim);
            } else {
                if (xo.hasAttribute(ROW_DIMENSION)) {
                    int rowDimension = xo.getIntegerAttribute(ROW_DIMENSION);
                    matrixParameter.setRowDimension(rowDimension);
                }

                if (xo.hasAttribute(COLUMN_DIMENSION)) {
                    int columnDimension = xo.getIntegerAttribute(COLUMN_DIMENSION);
                    matrixParameter.setColumnDimension(columnDimension);
                }
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
                AttributeRule.newIntegerRule(COLUMN_DIMENSION, true),
                AttributeRule.newBooleanRule(TRANSPOSE, true),
                AttributeRule.newBooleanRule(AS_COMPOUND, true),
                AttributeRule.newBooleanRule(BEHAVIOR, true),
        };

        public Class getReturnType() {
            return MatrixParameter.class;
        }
    };
}
