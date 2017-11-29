/*
 * FastMatrixParameter.java
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
import java.util.Arrays;
import java.util.List;

/**
 * @author Marc A. Suchard
 */
public class FastMatrixParameter extends CompoundParameter implements MatrixParameterInterface {

    public static final String FAST_MATRIX_PARAMETER = "fastMatrixParameter";
    public static final String ROW_DIMENSION = MatrixParameter.ROW_DIMENSION;
    public static final String COLUMN_DIMENSION = MatrixParameter.COLUMN_DIMENSION;
    public static final String STARTING_VALUE = "startingValue";

    public FastMatrixParameter(String id, int rowDimension, int colDimension, double startingValue) {
        super(id);
        singleParameter = new Parameter.Default(rowDimension * colDimension);
        addParameter(singleParameter);
        for (int i = 0; i < singleParameter.getDimension(); i++) {
            singleParameter.setParameterValue(i, startingValue);
        }
        this.rowDimension = rowDimension;
        this.colDimension = colDimension;
        DefaultBounds bounds = new DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, singleParameter.getDimension());
        addBounds(bounds);
    }

    public Parameter getParameter(int index) {
        if (proxyList == null) {
            proxyList = new ArrayList<ParameterProxy>(colDimension);
            for (int i = 0; i < colDimension; ++i) {
                proxyList.add(new ParameterProxy(this, i));
            }
        }
        return proxyList.get(index);
    }

    class ParameterProxy extends Parameter.Abstract {

        private final int column;
        private final FastMatrixParameter matrix;

        ParameterProxy(FastMatrixParameter matrix, int column) {
            this.matrix = matrix;
//            this.addParameterListener(this.matrix);
            this.column = column;
        }

        @Override
        protected void storeValues() {
            // Do nothing; storeValues() on whole matrix should have been called.
        }

        @Override
        protected void restoreValues() {
            // Do nothing; restoreValues() on whole matrix should have been called.
        }

        @Override
        protected void acceptValues() {
            // Do nothing; acceptValues() on whole matrix should have been called.
        }

        @Override
        protected void adoptValues(Parameter source) {
            throw new RuntimeException("Do not call");
        }

        @Override
        public double getParameterValue(int dim) {
            return matrix.getParameterValue(dim, column);
        }

        @Override
        public void setParameterValue(int dim, double value) {
            matrix.setParameterValue(dim, column, value);
            super.fireParameterChangedEvent(dim, Parameter.ChangeType.VALUE_CHANGED);
        }

        @Override
        public void setParameterValueQuietly(int dim, double value) {
            matrix.setParameterValueQuietly(dim, column, value);
        }

        @Override
        public void setParameterValueNotifyChangedAll(int dim, double value) {
            throw new RuntimeException("Do not call");
        }

        @Override
        public String getParameterName() {
            return getId();
        }

        @Override
        public void addBounds(Bounds<Double> boundary) {
            bounds = boundary;
        }

        @Override
        public Bounds<Double> getBounds() {
            return bounds;
        }

        @Override
        public void fireParameterChangedEvent(int index, ChangeType type){
            matrix.fireParameterChangedEvent(index + column * getDimension(), type);
            super.fireParameterChangedEvent(index, ChangeType.VALUE_CHANGED);
        }

        @Override
        public void addDimension(int index, double value) {
            throw new RuntimeException("Do not call");
        }

        @Override
        public double removeDimension(int index) {
            throw new RuntimeException("Do not call");
        }

        @Override
        public int getDimension() {
            return matrix.getRowDimension();
        }

        private Bounds<Double> bounds = null;

    }

    private final int index(int row, int col) {
        // column-major
        if(col > getColumnDimension()){
            throw new RuntimeException("Column " + col + " out of bounds: Compared to " + getColumnDimension() + "maximum size.");
        }
        if(row > getRowDimension()){
            throw new RuntimeException("Row " + row + " out of bounds: Compared to " + getRowDimension() + "maximum size.");
        }
        return col * rowDimension + row;
    }

    @Override
    public double getParameterValue(int row, int col) {
        return singleParameter.getParameterValue(index(row, col));
    }

    @Override
    public double[] getParameterValues() {
        double[] destination = new double[getDimension()];
        copyParameterValues(destination, 0);
        return destination;
    }

    @Override
    public void copyParameterValues(double[] destination, int offset) {
        final double[] source = ((Parameter.Default) singleParameter).inspectParameterValues();
        System.arraycopy(source, 0, destination, offset, source.length);
    }

    @Override
    public void setAllParameterValuesQuietly(double[] values, int offset) {
        final double[] destination = ((Parameter.Default) singleParameter).inspectParameterValues();
        System.arraycopy(values, offset, destination, 0, destination.length);
    }

    @Override
    public String toSymmetricString() {
        return MatrixParameter.toSymmetricString(this);
    }

    @Override
    public void setParameterValue(int row, int col, double value) {
        singleParameter.setParameterValue(index(row, col), value);
    }

    @Override
    public void setParameterValueQuietly(int row, int col, double value) {
        singleParameter.setParameterValueQuietly(index(row, col), value);
    }

    @Override
    public void setParameterValueNotifyChangedAll(int row, int col, double value) {
        singleParameter.setParameterValueNotifyChangedAll(index(row, col), value);
    }

    @Override
    public double[] getColumnValues(int col) {
        double[] rtn = new double[getRowDimension()];
        for (int i = 0; i < getRowDimension(); ++i) {
            rtn[i] = getParameterValue(i, col);
        }
        return rtn;
    }

    @Override
    public double[][] getParameterAsMatrix() {
        double[][] rtn = new double[getRowDimension()][getColumnDimension()];
        for (int j = 0; j < getColumnDimension(); ++j) {
            for (int i = 0; i < getRowDimension(); ++i) {
                rtn[i][j] = getParameterValue(i, j);
            }
        }
        return rtn;
    }

    @Override
    public int getColumnDimension() {
        return colDimension;
    }

    @Override
    public int getRowDimension() {
        return rowDimension;
    }

    @Override
    public int getParameterCount() {
        return getColumnDimension();
    }

    @Override
    public int getUniqueParameterCount() {
        return 1;
    }

    @Override
    public Parameter getUniqueParameter(int index) {
        return super.getParameter(0);
    }

    public void addBounds(Bounds<Double> boundary){
        singleParameter.addBounds(boundary);
    }

    public Bounds<Double> getBounds(){
        return singleParameter.getBounds();
    }


    private final int rowDimension;
    private final int colDimension;
    private final Parameter singleParameter;

    private List<ParameterProxy> proxyList = null;

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return FAST_MATRIX_PARAMETER;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            final String name = xo.hasId() ? xo.getId() : null;
            final double startingValue = xo.hasAttribute(STARTING_VALUE) ? xo.getDoubleAttribute(STARTING_VALUE) : 1;
            final int rowDimension = xo.getIntegerAttribute(ROW_DIMENSION);
            final int colDimension = xo.getIntegerAttribute(COLUMN_DIMENSION);

            FastMatrixParameter matrixParameter = new FastMatrixParameter(name, rowDimension, colDimension, startingValue);

            return matrixParameter;
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "A fast matrix parameter constructed from a single parameter.";
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                new ElementRule(Parameter.class, 0, Integer.MAX_VALUE),
                AttributeRule.newIntegerRule(ROW_DIMENSION, false),
                AttributeRule.newIntegerRule(COLUMN_DIMENSION, false),
                AttributeRule.newDoubleRule(STARTING_VALUE, true),
        };

        public Class getReturnType() {
            return FastMatrixParameter.class;
        }
    };
}
